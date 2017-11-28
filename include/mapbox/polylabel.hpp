#pragma once

#include <mapbox/geometry/polygon.hpp>
#include <mapbox/geometry/envelope.hpp>
#include <mapbox/geometry/point.hpp>
#include <mapbox/geometry/point_arithmetic.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>

namespace mapbox {

namespace detail {

// get squared distance from a point to a segment
template <class T>
T getSegDistSq(const geometry::point<T>& p,
               const geometry::point<T>& a,
               const geometry::point<T>& b) {
    auto x = a.x;
    auto y = a.y;
    auto dx = b.x - x;
    auto dy = b.y - y;

    if (dx != 0 || dy != 0) {

        auto t = ((p.x - x) * dx + (p.y - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = b.x;
            y = b.y;

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = p.x - x;
    dy = p.y - y;

    return dx * dx + dy * dy;
}

// signed distance from point to polygon outline (negative if point is outside)
template <class T>
auto pointToPolygonDist(const geometry::point<T>& point, const geometry::polygon<T>& polygon) {
    bool inside = false;
    auto minDistSq = std::numeric_limits<double>::infinity();

    for (const auto& ring : polygon) {
        for (std::size_t i = 0, len = ring.size(), j = len - 1; i < len; j = i++) {
            const auto& a = ring[i];
            const auto& b = ring[j];

            if ((a.y > point.y) != (b.y > point.y) &&
                (point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x)) inside = !inside;

            minDistSq = std::min(minDistSq, getSegDistSq(point, a, b));
        }
    }

    return (inside ? 1 : -1) * std::sqrt(minDistSq);
}

template <class T>
struct FitnessFunctor {
    FitnessFunctor(geometry::point<T> centroid, geometry::point<T> polygonSize)
        : centroid(centroid),
          maxSize(std::max(polygonSize.x, polygonSize.y))
        {}

    T operator()(const geometry::point<T>& cellCenter, T distancePolygon) const {
        if (distancePolygon <= 0) {
            return distancePolygon;
        }
        geometry::point<T> d = cellCenter - centroid;
        double distanceCentroid = std::sqrt(d.x * d.x + d.y * d.y);
        return distancePolygon * (1 - distanceCentroid / maxSize);
    }

    geometry::point<T> centroid;
    T maxSize;
};

template <class T>
struct Cell {
    template <class FitnessFunc>
    Cell(const geometry::point<T>& c_, T h_, const geometry::polygon<T>& polygon, const FitnessFunc& ff)
        : c(c_),
          h(h_),
          d(pointToPolygonDist(c, polygon)),
          fitness(ff(c, d)),
          maxFitness(ff(c, d + h * std::sqrt(2)))
        {}

    geometry::point<T> c; // cell center
    T h; // half the cell size
    T d; // distance from cell center to polygon
    T fitness; // fitness of the cell center
    T maxFitness; // a "potential" of the cell calculated from max distance to polygon within the cell
};

// get polygon centroid
template <class T>
geometry::point<T> getCentroid(const geometry::polygon<T>& polygon) {
    T area = 0;
    geometry::point<T> c { 0, 0 };
    const auto& ring = polygon.at(0);

    for (std::size_t i = 0, len = ring.size(), j = len - 1; i < len; j = i++) {
        const geometry::point<T>& a = ring[i];
        const geometry::point<T>& b = ring[j];
        auto f = a.x * b.y - b.x * a.y;
        c.x += (a.x + b.x) * f;
        c.y += (a.y + b.y) * f;
        area += f * 3;
    }

    return area == 0 ? ring.at(0) : c / area;
}

} // namespace detail

template <class T>
geometry::point<T> polylabel(const geometry::polygon<T>& polygon, T precision = 1, bool debug = false) {
    using namespace detail;

    // find the bounding box of the outer ring
    const geometry::box<T> envelope = geometry::envelope(polygon.at(0));

    const geometry::point<T> size {
        envelope.max.x - envelope.min.x,
        envelope.max.y - envelope.min.y
    };

    const T cellSize = std::min(size.x, size.y);
    T h = cellSize / 2;

    // a priority queue of cells in order of their "potential" (max distance to polygon)
    auto compareMax = [] (const Cell<T>& a, const Cell<T>& b) {
        return a.maxFitness < b.maxFitness;
    };
    using Queue = std::priority_queue<Cell<T>, std::vector<Cell<T>>, decltype(compareMax)>;
    Queue cellQueue(compareMax);

    if (cellSize == 0) {
        return envelope.min;
    }

    geometry::point<T> centroid = getCentroid(polygon);
    FitnessFunctor<T> fitnessFunc(centroid, size);

    // cover polygon with initial cells
    for (T x = envelope.min.x; x < envelope.max.x; x += cellSize) {
        for (T y = envelope.min.y; y < envelope.max.y; y += cellSize) {
            cellQueue.push(Cell<T>({x + h, y + h}, h, polygon, fitnessFunc));
        }
    }

    // take centroid as the first best guess
    auto bestCell = Cell<T>(centroid, 0, polygon, fitnessFunc);

    auto numProbes = cellQueue.size();
    while (!cellQueue.empty()) {
        // pick the most promising cell from the queue
        auto cell = cellQueue.top();
        cellQueue.pop();

        // update the best cell if we found a better one
        if (cell.fitness > bestCell.fitness) {
            bestCell = cell;
            if (debug) std::cout << "found best " << ::round(1e4 * cell.d) / 1e4 << " after " << numProbes << " probes" << std::endl;
        }

        // do not drill down further if there's no chance of a better solution
        if (cell.maxFitness - bestCell.fitness <= precision) continue;

        // split the cell into four cells
        h = cell.h / 2;
        cellQueue.push(Cell<T>({cell.c.x - h, cell.c.y - h}, h, polygon, fitnessFunc));
        cellQueue.push(Cell<T>({cell.c.x + h, cell.c.y - h}, h, polygon, fitnessFunc));
        cellQueue.push(Cell<T>({cell.c.x - h, cell.c.y + h}, h, polygon, fitnessFunc));
        cellQueue.push(Cell<T>({cell.c.x + h, cell.c.y + h}, h, polygon, fitnessFunc));
        numProbes += 4;
    }

    if (debug) {
        std::cout << "num probes: " << numProbes << std::endl;
        std::cout << "best distance: " << bestCell.d << std::endl;
    }

    return bestCell.c;
}

} // namespace mapbox
