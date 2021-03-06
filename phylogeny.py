# Phylogenetic trees in SAGE. When rooted, 0 is always the root, then 1 through
# n are the leaves and higher values are used for the internal nodes. When
# unrooted, 1 through n are the leaves. No promises are made concerning the
# numbering of the internal nodes.

# TODO: rooted is determined only by having a zero leaf?

import copy
from sage.all import Graph, matrix, Permutation, SymmetricGroup
from sage.plot.graphics import GraphicsArray


class Phylogeny(Graph):
    def __init__(self, rooted=True, *args, **kwargs):
        super(Phylogeny, self).__init__(*args, **kwargs)
        self.rooted = rooted

    # TODO: it's still not completely clear why I need to do this explicitly
    # for a subclass. Furthermore, there is this:
    # /home/matsen/re/curvature/tangle/all-hail-sage/phylogeny.py in copy(self, **kwargs)
    #      16     # for a subclass.
    #      17     def copy(self, **kwargs):
    # ---> 18         return copy.deepcopy(self, **kwargs)
    #      19
    #      20     def _repr_(self):
    #
    # TypeError: deepcopy() got an unexpected keyword argument 'immutable'
    #
    # so we absorb the kwargs here.
    def copy(self, **kwargs):
        return copy.deepcopy(self)

    def _repr_(self):
        return self.to_newick()

    def n_leaves(self):
        """
        Return the number of leaves.
        """
        # TODO: store in object?
        n = len(self.leaf_edges())
        if self.rooted:
            return n-1  # Root doesn't count.
        else:
            return n

    def plot(self):
        if self.rooted:
            return super(Phylogeny, self).plot(
                layout='tree', tree_root=0, tree_orientation="down")
        else:
            return super(Phylogeny, self).plot()

    def check_rooted_zero_edge(self):
        if not self.has_vertex(0):
            raise ValueError('Zero vertex not present in this rooted tree.')
        if self.order() > 1 and self.edges_incident(0) == []:
            raise ValueError('Lonely zero vertex!')

    def rooted_at(self, v):
        """
        Makes a new version of self rooted at v.
        """
        h = self.copy()
        if h.has_vertex(0):
            h.delete_vertex(0)
        h.add_edge(0, v)
        h.rooted = True
        return h

    def rooted_reduce(self, f_internal, f_leaf):
        """
        Assume that t is rooted at 0 and recur down through the rest of the
        tree using the supplied functions.
        """
        self.check_rooted_zero_edge()
        # Imagine an arrow pointing from src to dst. That is where the subtree
        # starts, with dst being the "root" node of the subtree.
        def aux(src, dst):
            nbs = self.neighbors(dst)
            nbs.remove(src)
            if nbs == []:  # Leaf.
                return f_leaf(dst)
            else:  # Internal node.
                return f_internal(
                        dst,
                        [aux(dst, daughter) for daughter in nbs])
        return aux(0, self.internal_root())

    def _rooted_to_newick(self, standalone=True):
        """
        Returns a Newick string such that the order of the subtrees is
        increasing in terms of the minimum of the leaf labels.
        `standalone` determines if it's part of a larger tree or not.
        """
        assert(self.rooted)
        # We carry along minimum leaf number as the first elt of a tuple to sort.
        def sorted_join(_, daughters):
            assert(daughters is not [])
            # Sort by this minimum leaf number.
            sorted_daughters = sorted(daughters, key=lambda d: d[0])
            return (
                sorted_daughters[0][0],  # The minimum leaf number.
                '('+','.join(d[1] for d in sorted_daughters)+')')
        self.check_rooted_zero_edge()
        if self.order() == 1:
            assert(standalone);  # Can't just have the root as a part of a tree.
            return '();'
        if self.order() == 2:
            # A 1-leaf rooted tree.
            if standalone:
                return '({});'.format(self.vertices()[1])
            return str(self.vertices()[1])
        _, nwk = self.rooted_reduce(sorted_join, lambda x: (x, str(x)))
        if standalone:
            return nwk+";"
        return nwk

    def to_newick(self):
        """
        Returns a Newick string such that the order of the subtrees is
        increasing in terms of the minimum of the leaf labels.
        """
        if self.rooted:
            return self._rooted_to_newick()
        else:
            # Root at the internal root.
            return self.rooted_at(self.internal_root())._rooted_to_newick()

    def _rooted_to_newick_shape(self, standalone=True):
        """
        Returns a Newick string such that the order of the subtrees is
        increasing in terms of lexicographical order.
        `standalone` determines if it's part of a larger tree or not.
        """
        assert(self.rooted)
        def sorted_join(_, daughters):
            return '('+','.join(sorted(daughters))+')'
        self.check_rooted_zero_edge()
        if self.order() == 1:
            assert(standalone);  # Can't just have the root as a part of a tree.
            return ';'
        if self.order() == 2:
            # A 1-leaf rooted tree.
            if standalone:
                return '();'
            return '()'
        nwk = self.rooted_reduce(sorted_join, lambda x: '')
        if standalone:
            return nwk+";"
        return nwk

    def to_newick_shape(self):
        """
        Return a Newick string representation of the shape (i.e. non-leaf-labeled
        graph rooted at zero) of tree t with larger trees always on the left.
        For unrooted trees, we root the tree at the centroid if it is unique.
        If not, then we take the centroid with the maximum-lexicographic shape
        string as defined by Python.
        """
        if self.rooted:
            return self._rooted_to_newick_shape()
        else:
            s = ""
            h = self.rooted_at(self.internal_root())
            for c in h.centroids():
                h.delete_edges(h.edges_incident(0))
                h.add_edge(0, c)
                new_s = h._rooted_to_newick_shape()
                if new_s > s:
                    s = new_s
            return s

    def make_special(self):
        """
        Make a graph out of the tree that is different than the original tree,
        and different from the result of make_special on any other distinct
        tree.
        """
        m = self.copy()
        m.allow_multiple_edges(True)
        if self.rooted:
            # Double up the edge to zero.
            m.add_edge(m.root_leaf(), m.internal_root())
        else:
            # Double up every leaf edge (no distinguished edges, so have to do
            # them all).
            for (u, v, _) in self.leaf_edges():
                m.add_edge(u, v)
        return m

    def rooted_is_isomorphic(self, other):
        """
        Are the two trees isomorphic if we give special status to the root?
        """
        return self.make_special().is_isomorphic(
            other.make_special(),
            certify=True)

    def leaf_edges(self):
        """
        Find all the edges corresponding to leaves in self.
        Here the "edge to the root" will count as a leaf.
        """
        return [(u, v, l) for (u, v, l) in self.edges()
                if self.degree(u) == 1 or self.degree(v) == 1]

    def root_leaf(self):
        """
        This is the root leaf for rooted trees, and an arbitrary choice of 1 for unrooted trees.
        """
        if self.rooted:
            return 0
        else:
            return 1

    def internal_root(self):
        """
        Return the node connecting to the zero leaf if rooted, or the one leaf if not.
        """
        assert(self.has_vertex(self.root_leaf()))
        n = self.neighbors(self.root_leaf())
        if n == []:
            raise ValueError('Empty tree has no internal root.')
        if len(n) > 1:
            raise ValueError('Too many neighbors of root leaf.')
        return n[0]

    def multiedge_leaf_edges(self):
        """
        Make a tree such that graph isomorphism is equivalent to leaf-labeled
        (labeled by the leaf vertex numbers) isomorphism.
        We assume that the tree was made by the enumeration code below, so that
        the leaf indices are all smaller than any internal node.
        """
        m = self.copy()
        m.allow_multiple_edges(True)
        for (u, v, _) in m.leaf_edges():
            # min(u, v) is the leaf number; we add the (leaf number + 1) edges
            # to a leaf edge so that it gets distinguished from other leaf
            # edges. Thus edge 0 will be doubled, edge 1 will be tripled, etc.
            for _ in range(min(u, v) + 1):
                m.add_edge(u, v)
        return m

    def llt_is_isomorphic(self, other, certificate=False):
        """
        Return if trees are leaf-labeled (labeled by the leaf vertex numbers)
        isomorphic.
        """
        return self.multiedge_leaf_edges().is_isomorphic(
            other.multiedge_leaf_edges(),
            certificate)

    def automorphism_group(self):
        """
        The automorphism group of the leaf nodes of a tree as a
        subgroup of the symmetric group.
        Like everything here, assumes that the first n nodes are leaves.
        """
        # If rooted, discount the root "leaf", and if unrooted this is the
        # correct max index:
        n = self.n_leaves()
        G = SymmetricGroup(n)
        if self.rooted:
            # Only take graph automorphisms that don't move 0.
            A = super(Phylogeny, self).automorphism_group(
                partition=[[0], range(1, self.order())])
        else:
            A = super(Phylogeny, self).automorphism_group()
        return G.subgroup(
            filter(
                # Just take the movement of the leaves, not the internals.
                lambda tup: all(i <= n for i in tup),
                g.cycle_tuples())
            for g in A.gens())

    def act_on_right(self, permish):
        """
        Returns a new tree.
        `permish` is something that has a .dict() method with the same number
        of values as there are leaves of the tree.
        """
        p = permish.dict().values()
        n = self.n_leaves()
        if n != len(p):
            raise ValueError(
                'Permutation must have same length as # leaves of tree.')
        # relabel has some "quirks":
        # http://www.sagemath.org/doc/reference/graphs/sage/graphs/generic_graph.html#sage.graphs.generic_graph.GenericGraph.relabel
        # Fix everything except for leaves.
        if self.rooted:
            p = [0]+p
        return self.relabel(p+range(n+1, self.order()), inplace=False)

    def centroids(self):
        """
        Returns centroids for trees with a zero leaf, ignoring that zero leaf.
        """
        # Tacky scoping workaround https://www.python.org/dev/peps/pep-3104/
        class Namespace:
            pass
        ns = Namespace()

        assert(self.has_vertex(0))
        n = len(self.leaf_edges()) - 1  # Always discount zero leaf.
        ns.best_score = n
        ns.candidates = [0]
        def aux_internal(curr, below_weights):
            tot = sum(below_weights)
            score = max([n - tot] + below_weights)
            if score < ns.best_score:
                ns.best_score = score
                ns.candidates = [curr]
            elif score == ns.best_score:
                ns.candidates.append(curr)
            return tot
        self.rooted_reduce(aux_internal, lambda _: 1)
        if n != 0:
            assert(ns.candidates != [0])  # Should find _something_ better!
        return ns.candidates


# Functions

def plot_tree_list(l):
    return GraphicsArray([t.plot() for t in l])


def _enumerate_rooted_bifurcating_trees(max_leaf_label, start_internal):
    """
    Rooted tree enumeration with internal nodes starting at some value and such
    that the maximum leaf label is `max_leaf_label`.
    """
    assert(max_leaf_label > 0)
    if max_leaf_label == 1:
        g = Phylogeny(rooted=True)
        g.add_vertices([0, 1])
        g.add_edge(0, 1)
        return [g]
    else:
        l = []
        # For every tree with one fewer taxon...
        for g in _enumerate_rooted_bifurcating_trees(
                                max_leaf_label-1, start_internal+1):
            # cut every edge in two and attach the new taxon, which is named
            # max_leaf_label.
            for (u, v, _) in g.edges():
                h = g.copy()
                h.delete_edge(u, v)
                h.add_vertex(max_leaf_label)
                # The new internal node is labeled with start_internal.
                h.add_vertex(start_internal)
                h.add_edges([(u, start_internal),
                             (v, start_internal),
                             (max_leaf_label, start_internal)])
                l.append(h)
        return l


def enumerate_bifurcating_trees(n_leaves, rooted=True):
    """
    Construct all the bifurcating, leaf-labeled phylogenetic trees.
    """
    if rooted:
        return _enumerate_rooted_bifurcating_trees(n_leaves, n_leaves+1)
    if n_leaves == 1:
        t = Phylogeny(rooted=rooted)
        return [t.add_vertices([1])]
    # Unrooted trees are in our setup are the same as rooted trees with one
    # less leaf (note the -1 below).
    l = _enumerate_rooted_bifurcating_trees(n_leaves-1, n_leaves)
    for t in l:
        # Boost labels by 1 so that leaves are 1 through n.
        t.relabel(range(1, t.order()+1))
        t.rooted = False
    return l


def indexed_tree_list(to):
    """
    Return a list of trees up to `to-1` such that the ith element in the list
    is a list of trees with i leaves.
    """
    return [[]] + [enumerate_rooted_trees(i) for i in range(1, to)]
