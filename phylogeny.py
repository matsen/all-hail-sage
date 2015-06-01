# Degree-3-internal-node phylogenetic trees in SAGE. When rooted, 0 is always
# the root, then 1 through n are the leaves and higher values are used for the
# internal nodes. When unrooted, 0 through n-1 are the leaves. No promises are
# made concerning the numbering of the internal nodes.

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

    def tree_reduce(self, f_internal, f_leaf):
        """
        Assume that t is rooted at 0 and recur down through the rest of the
        tree using the supplied functions.
        """
        # Imagine arrow pointing from src to dst.
        # That is where the subtree starts.
        def aux(src, dst):
            n = self.neighbors(dst)
            n.remove(src)
            if n == []:  # Leaf.
                return f_leaf(dst)
            else:  # Internal node.
                [left, right] = n
                return f_internal(aux(dst, left), aux(dst, right))
        return aux(0, self.internal_root())

    def _rooted_to_newick(self, standalone=True):
        """
        Returns a Newick string such that the order of the subtrees is
        increasing in terms of the minimum of the leaf labels.
        `standalone` determines if it's part of a larger tree or not.
        """
        # Carry along minimum leaf number to sort.
        def sorted_join((a, a_str), (b, b_str)):
            if a < b:
                return (a, '('+a_str+','+b_str+')')
            else:
                return (b, '('+b_str+','+a_str+')')
        if len(self.neighbors(0)) == 0:
            return '();'
        if self.order() == 2:
            # A 1-leaf rooted tree.
            if standalone:
                return '({});'.format(self.internal_root())
            return str(self.internal_root())
        _, nwk = self.tree_reduce(sorted_join, lambda x: (x, str(x)))
        if standalone:
            return nwk+";"
        return nwk

    def to_newick(self):
        if self.n_leaves() <= 2 or self.rooted:
            return self._rooted_to_newick()
        else:
            [z, l, r] = self.neighbor_trees(self.internal_root())
            # z is just the zero leaf, which we add back in via the string
            # below.
            assert(z.vertices() == [0])
            # TODO: sort subtrees.
            return '(0,'+\
                ','.join([
                    t._rooted_to_newick(standalone=False)
                    for t in [l, r]])+\
                ');'

    def _rooted_to_newick_shape(self, standalone=True):
        """
        Returns a Newick string such that the order of the subtrees is
        increasing in terms of lexicographical order.
        `standalone` determines if it's part of a larger tree or not.
        """
        def sorted_join(a, b):
            if a < b:
                return ('('+a+','+b+')')
            else:
                return ('('+b+','+a+')')
        if len(self.neighbors(0)) == 0 or self.order() == 2:
            # A 1-leaf rooted tree.
            if standalone:
                return '();'
            return ''
        nwk = self.tree_reduce(sorted_join, lambda _: "")
        if standalone:
            return nwk+";"
        return nwk

    def to_newick_shape(self):
        """
        Return a Newick string representation of the shape (i.e. non-leaf-labeled
        graph rooted at zero) of tree t with larger trees always on the left.
        """
        if self.n_leaves() <= 2 or self.rooted:
            return self._rooted_to_newick_shape()
        else:
            [z, l, r] = self.neighbor_trees(self.internal_root())
            # z is just the zero leaf, which we add back in via the string
            # below.
            assert(z.vertices() == [0])
            strings = sorted(map(
                            lambda t: t._rooted_to_newick_shape(standalone=False),
                            [l, r]))
            return '(,'+','.join(strings)+');'

    def duplicate_zero_edge(self):
        """
        Return a new tree that is the same as the original except with the edge
        to zero doubled up.
        """
        m = self.copy()
        m.allow_multiple_edges(True)
        [internal_root] = self.neighbors(0)
        m.add_edge(0, internal_root)
        return m

    def rooted_is_isomorphic(self, other):
        """
        Are the two trees isomorphic if we give special status to the root?
        """
        return self.duplicate_zero_edge().is_isomorphic(
            other.duplicate_zero_edge(),
            certify=True)

    def leaf_edges(self):
        """
        Find all the edges corresponding to leaves in self.
        Here the "edge to the root" will count as a leaf.
        """
        return [(u, v, l) for (u, v, l) in self.edges()
                if self.degree(u) == 1 or self.degree(v) == 1]

    def internal_root(self):
        """
        Return the node connecting to the zero leaf.
        """
        if not self.has_vertex(0):
            raise ValueError('This graph has no 0 vertex.')
        n = self.neighbors(0)
        if n == []:
            raise ValueError('Empty tree has no internal root.')
        if len(n) > 1:
            print n
            raise ValueError('0 vertex is not a leaf.')
        return n[0]

    def correct_deg_two_vertex(self, v):
        """
        If v is a degree two vertex, remove it and heal the edge.
        """
        n = self.neighbors(v)
        if len(n) == 2:
            self.delete_vertex(v)
            self.add_edge(n[0], n[1])

    def neighbor_trees(self, v):
        """
        Return the three trees around an internal node, each of which gets
        rooted at 0.
        """
        # TODO: further validation.
        g = self.copy()
        neighbors = self.neighbors(v)
        if len(neighbors) == 1:
            raise ValueError('`neighbor_trees` expects an internal node.')
        g.delete_vertex(v)
        def rooted_subtree_with_root(w):
            h = g.copy()
            h.subgraph(
                vertices=g.connected_component_containing_vertex(w), inplace=True)
            h.add_vertex(0)
            h.add_edge(0, w)
            h.correct_deg_two_vertex(w)
            h.rooted = True # TODO: eventually add a rooting method?
            return h
        return [rooted_subtree_with_root(w) for w in neighbors]

    def multiedge_leaf_edges(self):
        """
        Make a tree such that graph isomorphism is equivalent to leaf-labeled
        (labeled by the leaf vertex numbers) isomorphism.
        We assume that the tree was made by enumerate_rooted_trees, so that the
        leaf indices are all smaller than any internal node.
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
        max_idx = self.n_leaves()
        if self.rooted:
            # Only take graph automorphisms that don't move 0.
            G = SymmetricGroup(max_idx)
            A = super(Phylogeny, self).automorphism_group(
                partition=[[0], range(1, self.order())])
        else:
            # Use this rather than symmetric group so that we can have 0 get
            # moved around.
            G = PermutationGroup([[(0,1)],[tuple(range(max_idx+1))]])
            A = super(Phylogeny, self).automorphism_group()
        return G.subgroup(
            filter(
                # Just take the movement of the leaves, not the internals.
                lambda tup: all(i <= max_idx for i in tup),
                g.cycle_tuples())
            for g in A.gens())

    def act_on_right(self, permish):
        """
        Returns a new tree.
        `permish` is something that has a .dict() method.
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
            return self.relabel([0]+p+range(n+1, self.order()), inplace=False)
        else:
            return self.relabel(p+range(n, self.order()), inplace=False)


# Functions

def plot_tree_list(l):
    return GraphicsArray([t.plot() for t in l])


def _enumerate_rooted_trees(max_leaf_label, start_internal):
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
        for g in _enumerate_rooted_trees(max_leaf_label-1, start_internal+1):
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


def enumerate_trees(n_leaves, rooted=True):
    """
    Construct all the bifurcating, leaf-labeled phylogenetic trees.
    """
    if rooted == False:
        if n_leaves == 1:
            t = Phylogeny(rooted=False)
            return [t.add_vertices([0])]

        # Unrooted trees are in our setup are the same as rooted trees with one
        # less leaf (note the -1 below).
        l = _enumerate_rooted_trees(n_leaves-1, n_leaves)
        for t in l:
            t.rooted = False  # TODO: eventually add a rooting method?
        return l

    return _enumerate_rooted_trees(n_leaves, n_leaves+1)


def indexed_tree_list(to):
    """
    Return a list of trees up to `to-1` such that the ith element in the list
    is a list of trees with i leaves.
    """
    return [[]] + [enumerate_rooted_trees(i) for i in range(1, to)]
