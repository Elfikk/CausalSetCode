from numpy import random as nr
import numpy as np
import matplotlib.pyplot as plt
from Graph import Graph
import copy

def poisson_sprinkle(rho = 100, sort = True):
    # Twinkle twinkle little star. Specifically in 1D case, into the u,v 
    # (0,0) to (1,1) square. Performs a sprinkle at a density of rho, 
    # and returns the generated points ordered along u if sort is True.

    # Pick number of points to be generated using the Poisson distribution.
    N = nr.poisson(rho)

    # Sample N (u,v) points.
    points = nr.random_sample((N, 2))

    #Sort points based on value of u, as in null coordinates, if u_i < u_j
    #and v_i < v_j - sorting once makes this more efficient.
    if sort:
        points = np.array(sorted(points, key=lambda x: x[0]))

    return points

def poisson_sprinkle_rectangular(rho = 100, xlims = [0, 1], tlims = [0, 1],
                                 t_sort = True):
    #Sprinkle points into the provided rectangular bounds (whoop).

    # Pick number of points to be generated using the Poisson distribution.
    A = (tlims[1] - tlims[0]) * (xlims[1] - xlims[0])
    N = nr.poisson(rho * A) 

    # Sample N (x, t) points.
    points = nr.random_sample((N, 2))
    points[:,0] = points[:,0] * (xlims[1] - xlims[0]) + xlims[0]
    points[:,1] = points[:,1] * (tlims[1] - tlims[0]) + tlims[0]

    # Sort by the t-coordinate if wanted.
    if t_sort:
        points = np.array(sorted(points, key=lambda x: x[1]))

    return points

def find_relations(points):
    # Finds all order relations between points, given an array of (u,v)
    # points, sorted by u. 

    relations = {i: [] for i in range(len(points))}

    for i in range(len(points) - 1):
        for j in range(i + 1, len(points)):
            v_i = points[i][1]
            v_j = points[j][1]

            if v_i < v_j:
                relations[i].append(j)

    return relations

def find_relations_cylinder(points, L):
    #Finds all relations in Koluza-Klein cylinder, given (x,t) points
    #within a x in [0, L].

    relations = {i: [] for i in range(len(points))}

    for i in range(len(points) - 1):
        for j in range(i + 1, len(points)):
            x_i, t_i = points[i]
            x_j, t_j = points[j]

            t_c = t_i + 0.5 * L

            if t_j >= t_c:
                relations[i].append(j)
            else:
                delta_t = (t_j - t_i)**2

                # Corresponds to being in the original light cone.
                delta_x_base_case = (x_j - x_i)**2

                # Corresponds to being in the P' light cone.
                delta_x_case_2a = (x_j - x_i - L)**2

                # Corresponds to being in the P light cone.
                delta_x_case_2b = (x_j - x_i + L)**2

                if delta_t > delta_x_base_case or delta_t > delta_x_case_2a or \
                    delta_t > delta_x_case_2b:
                    relations[i].append(j)
    
    return relations

def find_relations_xt(points):
    #Finds all relations between points given x,t points.

    relations = {i: [] for i in range(len(points))}

    for i in range(len(points) - 1):
        for j in range(i + 1, len(points)):
            x_i, t_i = points[i]
            x_j, t_j = points[j]

            delta_t = (t_j - t_i)**2

            # Corresponds to being in the  light cone.
            delta_x= (x_j - x_i)**2

            if delta_t > delta_x:
                relations[i].append(j)
    
    return relations

def causal_links(points, relations = None):
    # Finds all order links between points, given an array of (u,v)
    # points, sorted by u. 
    # Is this really more efficient, or this just a nice trick?

    if not relations:
        relations = find_relations(points)

    links = {i: [] for i in range(len(points))}

    # print(relations)

    for i in range(len(relations)):
        precedes = set(relations[i].copy())
        to_remove = set()
        for j in precedes:
            to_remove.update(relations[j])
        links_set = precedes.difference(to_remove)
        links[i] = list(links_set)

    return links

def spacelike_relations(points):
    # Expects points in ascending order of u.

    # points = points[::-1]

    spacelike_relations = {i: [] for i in range(len(points))}

    for i in range(len(points) - 1):
        for j in range(i + 1, len(points)):
            v_i = points[i][1]
            v_j = points[j][1]

            if v_i > v_j:
                spacelike_relations[i].append(j)

    return spacelike_relations

def spacelike_links(points):
    # Finds space-like "links" - the smallest space-like relations between 
    # space-like points, to be used for manageable sized graphs for finding
    # longest maximal chain. 

    relations = spacelike_relations(points)

    links = {i: [] for i in range(len(points))}

    for i in range(len(relations)):
        precedes = set(relations[i].copy())
        to_remove = set()
        for j in precedes:
            to_remove.update(relations[j])
        links_set = precedes.difference(to_remove)
        links[i] = list(links_set)

    return links

def plot_links(points, links, colour = "black", style = "o--"):
    # Plots all links, given a list of points ordered along u and
    # the links of the causal set.
    
    plt.plot(points[:, 0], points[:, 1], "x")

    for i in range(len(links)):
        for j in links[i]:
            u = [points[i][0], points[j][0]]
            v = [points[i][1], points[j][1]]
            plt.plot(u, v, style, c = colour)

    plt.grid()
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)
    # plt.show()

def plot_links_xt(x, t, links, colour = "black", style = "o--"):

    # plt.plot(t, x, "x")

    for i in range(len(links)):
        for j in links[i]:
            xs = [x[i], x[j]]
            ts = [t[i], t[j]]
            plt.plot(xs, ts, style, c = colour, alpha = 0.4)

    # for i in range(len(x)):
    #     plt.plot([x[i], x[i] + 0.1], [t[i], t[i] + 0.1], "bo--", c = "blue")
    #     plt.plot([x[i] - 0.1, x[i]], [t[i] + 0.1, t[i]], "bo--", c = "blue")

    # for k in range(-1, 2):
    #     plt.plot([x[0] + k, x[0] + 0.5 + k], [t[0], t[0] + 0.5], "bo--")
    #     plt.plot([x[0] - 0.5 + k, x[0] + k], [t[0] + 0.5, t[0]], "bo--")
    

    # #Light cone
    plt.plot((0, 1), (0, 1), "bo--")
    plt.plot((0, -1), (0, 1), "bo--")

    plt.grid()
    # plt.xlim(-0.75, 0.75)
    # plt.ylim(-0.1, 1.5)

    plt.xlabel("x")
    plt.ylabel("t")

    plt.plot()
    
def plot_links_cylinder(x, t, links, L = 1, T = 10, colour = "black",
                        style = "o--"):
    

    for i in range(len(links)):
        for j in links[i]:
            xs = [x[i], x[j]]
            ts = [t[i], t[j]]
            plt.plot(xs, ts, style, c = colour, alpha = 0.4, markersize = 0)

    plt.plot(x, t, "x", c = "blue")

    # Cylinder Boundaries
    plt.plot((0, 0), (0, T), "bo--")
    plt.plot((L, L), (0, T), "bo--")

    plt.grid()
    plt.xlim(-0.1 * L, 1.1 * L)
    plt.ylim(-0.1 * T, 1.1 * T)

    plt.xlabel("x")
    plt.ylabel("t")

    # plt.plot()

def periodic_cylinder_link(xi, ti, xj, tj, L = 1, colour = "black",
                        style = "o--"):

    # Improvement - categorise relation/link type when checking for relation/link
    # in the first place rather than working it out here.

    time_comp = (tj - ti)**2

    x_comp_base = (xj - xi)**2
    
    
    if time_comp >= x_comp_base:
        xs = [xi, xj]
        ts = [ti, tj]
        plt.plot(xs, ts, style, c = colour, alpha = 0.4, markersize = 0)
    else:
        x_comp_2b = (xj - xi + L)**2
        if time_comp >= x_comp_2b:
            t_trans = ((tj - ti) / (xj + L - xi)) * (L - xi) + ti
            x1s = [xi, L]
            t1s = [ti, t_trans]
            x2s = [0, xj]
            t2s = [t_trans, tj]
        else:
            t_trans = ((tj - ti) / (xi + L - xj)) * xi + ti
            x1s = [xi, 0]
            t1s = [ti, t_trans]
            x2s = [L, xj]
            t2s = [t_trans, tj]

            plt.plot(x1s, t1s, style, c = colour, alpha = 0.4, markersize = 0)
            plt.plot(x2s, t2s, style, c = colour, alpha = 0.4, markersize = 0)


def plot_links_cylinder_periodic(x, t, links, L = 1, T = 10, 
                                 colour = "black", style = "o--"):

    for i in range(len(links)):
        for j in links[i]:
            xi, xj = x[i], x[j]
            ti, tj = t[i], t[j]
            periodic_cylinder_link(xi, ti, xj, tj, L, colour, style)

    plt.plot(x, t, "x", c = "blue")

    # Cylinder Boundaries
    plt.plot((0, 0), (0, T), "bo--",  markersize = 0)
    plt.plot((L, L), (0, T), "bo--",  markersize = 0)

    plt.grid()
    plt.xlim(-0.1 * L, 1.1 * L)
    plt.ylim(-0.1 * T, 1.1 * T)

    plt.xlabel("x")
    plt.ylabel("t")


def to_xt(u, v):
    #Converts u, v points to x, t.

    x = np.sqrt(2)/2 * (u + v)
    t = np.sqrt(2)/2 * (v - u)

    return x, t

def to_uv(x, t):
    # Converts x,t points to u,v.

    u = np.sqrt(2)/2 * (t - x)
    v = np.sqrt(2)/2 * (x + t)

    return u, v

def max_boundary(xi, x_c, t_c, h):
    if xi < x_c:
        t_max = xi - x_c + h + t_c
        return t_max
    else:
        t_max = x_c + h - xi + t_c
        return t_max

def min_boundary(xi, x_c, t_c, h):
    if xi < x_c:
        t_min = x_c - h - xi + t_c
        return t_min
    else:
        t_min = xi - x_c - h + t_c
        return t_min

def aleksandrov_interval_sample(points, l, L, T):
    # Given the list of all points within a cylindrical interval, takes a 
    # aleksandrov interval sample of a given side length l, with a random
    # centre point x_c, t_c. Returns points in aleksandrov interval.

    point_subset = set()

    # Half the height of the aleksandrov along constant t/x.
    h = l / np.sqrt(2)

    # Pick x coordinate. Periodic boundary so can pick any x point and have
    # interval overlap.
    x_c = L * np.random.random()

    # Pick t coordinate. Needs to not have the top or bottom overflowing.
    # So can range between h, T - h.
    t_c = (T - 2 * h) * np.random.random() + h

    # print(x_c, t_c)

    # The left-most and right-most x-coordinates of the aleksandrov interval
    x_l, x_r = x_c - h, x_c + h
    
    # print(x_l, x_r)

    # Copy points to work on directly but not affect points array outside of 
    # the function.
    points_copy = np.copy(points)

    #If the entire interval overlaps 0 to L, then don't shift anything.
    if x_l < 0 and x_r > L:
        pass
    #If x_l < 0, shift all points within L - x_l to L by -L
    elif x_l < 0:
        for p in points_copy:
            x_p = p[0]
            if x_p > L + x_l:
                p[0] -= L 
    #If x_r > L, shift all points within 0 to x_r - L by +L 
    elif x_r > L:
        for p in points_copy:
            x_p = p[0]
            if x_p < x_r - L:
                p[0] += L

    # Now apply standard checking for the case of interval somewhere within 
    # boundaries.

    for i in range(len(points)):
        
        xi = points_copy[i][0]
        #If within the aleksandrov interval, s.t x_l < xi < x_r (not equals
        #as that would have to fall in on the tip of the interval).
        if xi > x_l and xi < x_r:
            # Check if worth computing the min and max allowed times for given
            # xi. So if t_c - h < t_i < t_c + h, continue.
            ti = points[i][1]
            if t_c - h < ti and ti < t_c + h:
                #tc falls in a good enough range. 
                t_max = max_boundary(xi, x_c, t_c, h)
                t_min = min_boundary(xi, x_c, t_c, h)

                #If within interval, add point to the subset.
                if t_min < ti and ti < t_max:
                    # print(t_min, t_c, t_max)
                    point_subset.add(i)

    return point_subset, x_c, t_c

def generate_position_subset(points, subset):
    # Takes the subset of point IDs and returns the actual points for quick 
    # reference.

    position_subset = np.array([points[i] for i in subset])
    return position_subset    

def generate_relations_subset(relations, points_subset):
    # Takes a subset of points and filters through existing relations, only
    # returning the subset of relations of included points.

    # Construct a new dictionary of relations with the correct starting
    # points included.
    relations_subset = {num: [] for num in points_subset}

    # Pass through relations subset lists.
    for i in relations_subset:
        rel_copy = []

        for j in relations[i]:
            if j in points_subset:
                rel_copy.append(j)
        relations_subset[i] = rel_copy
    
    return relations_subset

if __name__ == "__main__":
    # p = poisson_sprinkle(100)
    # # print(len(p))
    # sp_links = spacelike_links(p)

    # # print(sp_links)
    # N = 0

    # for i in range(len(sp_links)):
    #     N += len(sp_links[i])

    # # print(N)

    # links = causal_links(p)    
    # print(links)

    # N = 0

    # for i in range(len(links)):
    #     N += len(links[i])

    # # print(N)

    # # plot_links(p, links)
    # plot_links(p, sp_links,colour = "red", style = "o-.")

    # u, v = p[:, 0], p[:, 1]

    # x, t = to_xt(u, v)

    # # plot_links_xt(x, t, sp_links, colour = "red", style = "o-.")
    # # plot_links_xt(x, t, links)
    # # plt.show()

    # spacelike_graph = Graph(len(sp_links))

    # for i in sp_links:
    #     for j in sp_links[i]:
    #         spacelike_graph.addEdge(i, j)

    # print(spacelike_graph.LongestPathLength())

    # # plot_links_xt(x, t, sp_links, colour = "red", style = "o-.")
    # # plot_links_xt(x, t, links)
    # plt.show()

    L = 1
    T = 2
    N_bar = 20
    rho = N_bar / (L * T)

    p = poisson_sprinkle_rectangular(rho, [0,L], [0,T])

    # print(p)

    relations = find_relations_cylinder(p, L)
    # print(relations[0])

    links = causal_links(p, relations)
    # print(links)

    x, t = p[:, 0], p[:, 1]

    plt.figure(figsize=(L*10, T*2.5))
    plot_links_cylinder_periodic(x, t, links, L, T)
    plt.show()