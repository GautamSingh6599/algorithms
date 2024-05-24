#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <ranges>
#include <stack>
#include <utility>
#include <vector>
// Copyright 2004 Gautam Singh

struct Point {
  int x, y;
};

// Find the Euclidean Distance between two points u, v.
double euclidean_distance(Point u, Point v) {
  return std::sqrt(std::pow(u.x - v.x, 2) + std::pow(u.y - v.y, 2));
}

double perimeter(Point polygon[], int i, int j, int k) {
  Point u = polygon[i], v = polygon[j], w = polygon[k];
  return euclidean_distance(u, v) + euclidean_distance(v, w) +
         euclidean_distance(w, u);
}

std::vector<std::vector<int>> minimum_cost_polygon_dp(Point polygon[], int n) {
  if (n < 3) {
    return {{0}};
  }
  std::vector<std::vector<double>> memo_t(n, std::vector<double>(n, -1));
  std::vector<std::vector<int>> k_store(n, std::vector<int>(n, -1));
  for (int gap = 0; gap < n; gap++) {
    for (int i = 0, j = gap; j < n; ++i, ++j) {
      if (j < i + 2) {
        memo_t[i][j] = 0.0;
      } else {
        memo_t[i][j] = std::numeric_limits<double>::max();
        for (int k = i + 1; k < j; ++k) {
          double cost =
              memo_t[i][k] + memo_t[k][j] + perimeter(polygon, i, j, k);
          if (memo_t[i][j] > cost) {
            memo_t[i][j] = cost;
            k_store[i][j] = k;
          }
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      std::cout << " |" << memo_t[i][j] << "| ";
    }
    std::cout << "\n";
  }

  std::cout << "------------" << "\n";
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      std::cout << " |" << k_store[i][j] << "| ";
    }
    std::cout << "\n";
  }

  std::cout << memo_t[0][n - 1] << "\n";
  return k_store;
}

std::vector<std::vector<int>>
Print_Triangulation(const std::vector<std::vector<int>> &k_store, int i,
                    int j) {
  std::vector<std::vector<int>> triangle;
  std::stack<std::pair<int, int>> s;
  s.push({i, j});
  while (!s.empty()) {
    auto [i, j] = s.top();
    s.pop();
    if (std::abs(i - j) > 1) {
      int k = k_store[i][j];
      triangle.push_back({i, k, j});
      s.push({k, j});
      s.push({i, k});
    }
  }
  return triangle;
}

double spanning_ratio() {
  double span_ratio;
  return span_ratio;
}

void addEdge(std::vector<std::vector<std::pair<int, double>>> &adj, int i,
             int j, double weight) {
  adj[i].push_back({j, weight});
  adj[j].push_back({i, weight});
}

int main() {
  Point polygon[] = {{0, 0}, {1, 0}, {2, 1}, {1, 2}, {0, 2}};
  int n = sizeof(polygon) / sizeof(polygon[0]);
  auto triangulation = minimum_cost_polygon_dp(polygon, n);
  auto triangles = Print_Triangulation(triangulation, 0, n - 1);
  for (const auto &i : triangles) {
    std::cout << "Triangles: ";
    for (int j : i) {
      std::cout << j << " ";
    }
    std::cout << "\n";
  }

  std::vector<std::vector<std::pair<int, double>>> adj[n];
  for (int i = 0; i < n - 1; i++) {
    addEdge(adj, i, i + 1, euclidean_distance(polygon[i], polygon[i + 1]));
  }

  for (auto &i : triangles) {
    std::ranges::sort(i);
    addEdge(adj, i[0], i[2], euclidean_distance(polygon[i[0]], polygon[i[2]]));
  }
  std::cout << "---------" << "\n";
  for (int i = 0; i < n; i++) {
    std::cout << "{ " << polygon[i].x << ", " << polygon[i].y << "} " << "\n";
    for (const auto &j : adj[i]) {
      for (auto [k, l] : j) {
        std::cout << "\t" << "{ " << polygon[k].x << ", " << polygon[k].y
                  << "} -----> " << l << "\n";
      }
    }
  }

  std::cout << spanning_ratio();
  return 0;
}
