#include <algorithm>
#include <cerrno>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>
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

std::vector<std::vector<std::pair<double, int>>>
minimum_cost_polygon_dp(Point polygon[], int n) {
  std::vector<std::vector<std::pair<double, int>>> memo(
      n, std::vector<std::pair<double, int>>(n));
  // if (n < 3) {
  //   return;
  // }
  for (int i = 0; i < n; i++) {
    for (int j = 0, k = i; k < n; j++, k++) {
      if (k < j + 2) {
        memo[j][k].first = 0.0;
        memo[j][k].second = -1;
      } else {
        memo[j][k].first = std::numeric_limits<double>::max();
        for (int l = j + 1; l < k; l++) {
          double cost =
              memo[j][l].first + memo[l][k].first + perimeter(polygon, j, k, l);
          if (memo[j][k].first > cost) {
            memo[j][k].first = cost;
            memo[j][k].second = l;
          }
        }
      }
    }
  }
  for (const auto &i : memo) {
    for (const auto &j : i) {
      std::cout << " |" << j.first << ", " << j.second << "| ";
    }
    std::cout << "\n";
  }
  return memo;
}

void Print_triangulation(
    const std::vector<std::vector<std::pair<double, int>>> &table, int i,
    int j) {
  if (i != j) {
    std::cout << i - 1 << " " << table[i][j].second << " " << j << "\n";
    Print_triangulation(table, i, table[i][j].second);
    Print_triangulation(table, table[i][j].second + 1, j);
  }
}

int main() {
  Point polygon[] = {{0, 0}, {1, 0}, {2, 1}, {1, 2}, {0, 2}};
  int n = sizeof(polygon) / sizeof(polygon[0]);
  auto triangulation = minimum_cost_polygon_dp(polygon, n);
  Print_triangulation(triangulation, 0, n - 1);
  return 0;
}
