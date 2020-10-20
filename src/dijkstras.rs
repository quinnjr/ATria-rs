// Copyright (C) Joseph R. Quinn
// SPDX-License-Identifier: MIT

pub fn dijkstras(graph: &Vec<Vec<f32>>, source: usize) -> Vec<f32> {
    let size = graph.len();
    let mut distances = vec![f32::INFINITY; size];
    let mut is_processed = vec![false; size];

    distances[source] = 0.0;

    for _count in 0..size - 1 {

        let u = min_distance(&distances, &is_processed);

        is_processed[u] = true;

        for v in 0..size {
            if !is_processed[v] &&
                graph[u][v] != f32::INFINITY &&
                distances[u] + graph[u][v] < distances[v] {
                distances[v] = distances[u] + graph[u][v];
            }
        }
    }

    distances
}

#[inline]
fn min_distance(distance: &Vec<f32>, processed: &[bool]) -> usize {
    let mut min = f32::INFINITY;
    let mut index = 0;

    for i in 0..distance.len() {
        if !processed[i] && distance[i] <= min {
            min = distance[i];
            index = i;
        }
    }

    index
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_werks() {
        let graph = vec![
            vec![0.0, f32::INFINITY, -2.0, f32::INFINITY],
            vec![4.0, 0.0, 3.0, f32::INFINITY],
            vec![f32::INFINITY, f32::INFINITY, 0.0, 2.0],
            vec![f32::INFINITY, -1.0, f32::INFINITY, 0.0]
        ];

        let expected = vec![
            vec![0.0, -1.0, -2.0, 0.0],
            vec![4.0, 0.0, 3.0, 5.0],
            vec![5.0, 1.0, 0.0, 2.0],
            vec![3.0, -1.0, 2.0, 0.0]
        ];

        let mut calculated = vec![];

        for i in 0..graph.len() {
            calculated.push(dijkstras(&graph, i));
        }

        assert_eq!(expected, calculated);
    }
}
