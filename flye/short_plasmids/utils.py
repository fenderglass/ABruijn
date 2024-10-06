#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

def find_connected_components(graph):
    def dfs(start_vertex, connected_components_counter):
        dfs_stack = [start_vertex]
        used[start_vertex] = True
        while len(dfs_stack):
            vertex = dfs_stack.pop()
            connected_components[vertex] = connected_components_counter
            for neighbour in graph[vertex]:
                if not used[neighbour]:
                    dfs_stack.append(neighbour)
                    used[neighbour] = True

    n_vertices = len(graph)
    connected_components = [0 for _ in range(n_vertices)]
    connected_components_counter = 0
    used = [False for _ in range(n_vertices)]

    for i in range(n_vertices):
        if not used[i]:
            dfs(i, connected_components_counter)
            connected_components_counter += 1

    return connected_components, connected_components_counter
