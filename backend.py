#!/usr/bin/env python2
# -*- coding: utf8 -*-

from fractions import gcd
from itertools import product, combinations, groupby

def enumerate_bare(Ne, M, d, Ct):
    """enumerate bare configurations

    Ne = nubmer of electrons
    M = number of spinless orbitals
    d = minimal distance between electrons. form cluster when distance = d
    Ct = maximal number of electrons in a given cluster

    return = a list of configs, where each config is a tuple of clusters,
             and each cluster is (j, size), where
                j  = the starting point of the cluster
                size = the number of electrons in the cluster
    """
    def enum(config):
        n_sofar = sum(size for j, size in config)
        if n_sofar == Ne:
            return [config]

        jfirst = M if len(config) == 0 else config[0][0]
        jnext = 0 if len(config) == 0 else (config[-1][0] + (config[-1][1] - 1) * d) + d + 1

        clusters = [(j, size)
                    for j in xrange(jnext, M)
                    for size in xrange(1, min(Ct, Ne - n_sofar) + 1)
                    if (j + (size - 1) * d) + d + 1 - M <= jfirst]
        return sum((enum(config + (cluster,)) for cluster in clusters), [])
    return enum(tuple())

def enumerate_dressed(configs, svalues):
    """enumerate dressed states

    configs = output of enumerate_bare
    svalues = [[s for i in xrange(Ct)] for j in xrange(M)],
              svalues[j] is the sorted list of the allowed values of s for a given j

    return = list of dressed states,
             where each dressed state is a dict {cluster: sorted tuple of s values}
    """
    return [dict(zip(config, ssets))
            for config in configs
            for ssets in product(*[combinations(svalues[j], size) for j, size in config])]

def count_orbits(states, Nx, Ny, C):
    """count the orbit lengths under the color-entangled magnetic translations

    states = a ky-subgroup of the output of enumerate_dressed
    NyCt = Ny / Ct

    return = [(length of orbit, sign), ...]
    """
    Ct = gcd(C, Ny)
    d = C / Ct
    M = (Nx * Ny) / Ct

    def T(state):
        """apply T: |j, s>  -->  |j + Ny/Ct, s + 1> to a dressed state

        return = (result of T on |j, s>, sign from antisymmetrization)
        """
        sign = 0
        result = {}
        for (j, size), sset in state.iteritems():
            jT = j + Ny / Ct

            ssetT = [(s + 1) - Nx * ((jT - (jT % M)) / M) for s in sset]

            ds = (ssetT[0] % C) - ssetT[0]
            n = sum(s + ds >= C for s in ssetT)
            sign += (size - n) * n

            result[(jT % M, size)] = tuple(sorted(s % C for s in ssetT))

        return result, (sign % 2)

    touched = [False] * len(states)
    orbits = []
    for i, initial in enumerate(states):
        if touched[i] is False:
            touched[i] = True
            size = 1
            state, sign = T(initial)

            while state != initial:
                size += 1

                index = states.index(state)
                assert touched[index] is False
                touched[index] = True

                state, sgn = T(state)
                sign += sgn

            sign %= 2
            if sign:
                assert Nx % 2 == 0

            orbits.append((size, sign))

    return orbits

def count(Ne, Nx, Ny, C):
    """count the number of zero modes for Ne bosons on Nx by Ny lattice with Chern number C

    return = {(kx, ky): degeneracy}
    """

    if Ne * (C + 1) > Nx * Ny:
        return dict()

    Ct = gcd(C, Ny)
    d = C / Ct
    M = (Nx * Ny) / Ct

    def shift_back(j, s):
        s = (s - Nx * ((j - (j % M)) / M)) % C
        j = j % M
        return j, s
    def get_js(X, ky):
        j = (X * Ny + ky * C) / Ct
        s = X % C
        return shift_back(j, s)

    js_list = [get_js(X, ky) for X in xrange(Nx) for ky in xrange(Ny)]
    assert len(js_list) == len(set(js_list))

    svalues = [sorted([s for j, s in js_list if j == j0]) for j0 in xrange(M)]

    configs = enumerate_bare(Ne, M, d, Ct)
    states = enumerate_dressed(configs, svalues)

    get_ky = lambda j, s: js_list.index(shift_back(j, s)) % Ny
    total_ky = lambda state: sum(get_ky(sum(xrange(j, j + d * size, d)), sum(sset))
                                 for (j, size), sset in state.iteritems()) % Ny
    states.sort(key=total_ky)

    counting = dict()
    for ky, group in groupby(states, total_ky):
        for n, sign in count_orbits(list(group), Nx, Ny, C):
            for kx in xrange(0, Nx, Nx / n):
                if sign:
                    kx = (kx + Nx % 2) % Nx
                counting.setdefault((kx, ky), 0)
                counting[(kx, ky)] += 1
    return counting

def pretty_print(counting, Nx, Ny, xlabel='x', ylabel='y', show_axes=True):
    """pretty print counting

    counting = {(kx, ky): c, }
    Nx, Ny = the size of BZ rectangle
    """
    output = []

    # field widths
    l = len(str(Ny - 1)) # width of ky ticks
    w = len(str(max(Nx - 1, max(counting.values() + [0])))) # width of each entry

    if show_axes:
        output.append(('k' + ylabel).ljust(l + 2) + ('tot=%s' % str(sum(counting.values()))).center((w + 2) * Nx + len('k' + xlabel)))
    for ky in reversed(xrange(Ny)):
        line = [str(ky).ljust(l + 1), u'├'] if show_axes else [] # ky axis
        for kx in xrange(Nx):
            if (kx, ky) in counting:
                c = counting[(kx, ky)]
                line.append(str(c).center(w + 2))
            else:
                line.append(u'·'.center(w + 2))
        output.append(''.join(line) + (' ' * len('k' + xlabel)) * show_axes)

    if show_axes:
        # kx axis
        line = (' ' * (l + 1)) + u'╰' + (u'┴'.center(w + 2, u'─') * Nx) + (' ' * len('k' + xlabel))
        output.append(line)

        # kx axis label
        line = [' ' * (l + 2)]
        for kx in xrange(Nx):
            line.append(str(kx).center(w+2))
        line.append('k' + xlabel)
        output.append(''.join(line))

    # need to encode to utf-8 explicitly for the special characters to work in a shell pipe
    # http://stackoverflow.com/a/492711
    for line in output:
        print line.encode('utf-8')
