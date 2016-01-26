import os, sys
import numpy as np
from scipy.io import mmread,mmwrite
from scipy.sparse import coo_matrix

def _transition_counts(sequences, n_states = 0, lag_time=1, sliding_window=True):
    """Count the number of directed transitions in a collection of sequences
    in a discrete space.
    Parameters
    ----------
    sequences : list of array-like
        List of sequences, or a single sequence. Each sequence should be a
        1D iterable of state labels. Labels can be integers, strings, or
        other orderable objects.
    lag_time : int
        The time (index) delay for the counts.
    sliding_window : bool
        When lag_time > 1, consider *all*
        ``N = lag_time`` strided sequences starting from index
         0, 1, 2, ..., ``lag_time - 1``. The total, raw counts will
         be divided by ``N``. When this is False, only start from index 0.
    Returns
    -------
    counts : array, shape=(n_states, n_states)
        ``counts[i][j]`` counts the number of times a sequences was in state
        `i` at time t, and state `j` at time `t+self.lag_time`, over the
        full set of trajectories.
    mapping : dict
        Mapping from the items in the sequences to the indices in
        ``(0, n_states-1)`` used for the count matrix.
    Examples
    --------
    >>> sequence = [0, 0, 0, 1, 1]
    >>> counts, mapping = _transition_counts([sequence])
    >>> print counts
    [[2, 1],
     [0, 1]]
    >>> print mapping
    {0: 0, 1: 1}
    >>> sequence = [100, 200, 300]
    >>> counts, mapping = _transition_counts([sequence])
    >>> print counts
    [[ 0.  1.  0.]
     [ 0.  0.  1.]
     [ 0.  0.  0.]]
    >>> print mapping
    {100: 0, 200: 1, 300: 2}
    Notes
    -----
    `None` and `NaN` are recognized immediately as invalid labels. Therefore,
    transition counts from or to a sequence item which is NaN or None will not
    be counted. The mapping return value will not include the NaN or None.
    """
    if (not sliding_window) and lag_time > 1:
        return _transition_counts([X[::lag_time] for X in sequences], lag_time=1)

    classes = np.unique(np.concatenate(sequences))
    contains_nan = (classes.dtype.kind == 'f') and np.any(np.isnan(classes))
    contains_none = any(c is None for c in classes)

    if contains_nan:
        classes = classes[~np.isnan(classes)]
    if contains_none:
        classes = [c for c in classes if c is not None]
    if n_states == 0:
        n_states = classes.max()+1

    mapping = dict(zip(classes, range(n_states)))
    mapping_is_identity = np.all(classes == np.arange(n_states))
    mapping_fn = np.vectorize(mapping.get, otypes=[np.int])
    none_to_nan = np.vectorize(lambda x: np.nan if x is None else x,
                               otypes=[np.float])

    counts = np.zeros((n_states, n_states), dtype=float)
    _transitions = []

    for y in sequences:
        y = np.asarray(y)
        from_states = y[: -lag_time: 1]
        to_states = y[lag_time::1]

        if contains_none:
            from_states = none_to_nan(from_states)
            to_states = none_to_nan(to_states)

        if contains_nan or contains_none:
            # mask out nan in either from_states or to_states
            mask = ~(np.isnan(from_states) + np.isnan(to_states))
            from_states = from_states[mask]
            to_states = to_states[mask]

        #if (not mapping_is_identity) and len(from_states) > 0 and len(to_states) > 0:
        from_states = mapping_fn(from_states)
        to_states = mapping_fn(to_states)

        _transitions.append(np.row_stack((from_states, to_states)))

    transitions = np.hstack(_transitions)
    C = coo_matrix((np.ones(transitions.shape[1], dtype=int), transitions),
        shape=(n_states, n_states))
    counts = counts + np.asarray(C.todense())

    # If sliding window is False, this function will be called recursively
    # with strided trajectories and lag_time = 1, which gives the desired
    # number of counts. If sliding window is True, the counts are divided
    # by the "number of windows" (i.e. the lag_time). Count magnitudes
    # will be comparable between sliding-window and non-sliding-window cases.
    # If lag_time = 1, sliding_window makes no difference.
    counts /= float(lag_time)

    return counts, mapping

if __name__ == "__main__":
    for i in range(6383,6391):
        fn = "../test/Assignments/Assignments-%d.npy"%i
        a = np.load(fn)
        counts,mapping = _transition_counts(a, n_states = 1000, lag_time=50)
        counts = coo_matrix(counts)
        mmwrite("tCounts-micro-%d.mtx"%i,counts)

