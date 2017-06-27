import os
import sys
import heapq

try:
    from pybloomfilter import BloomFilter
except ImportError:
    raise ImportError('pybloomfilter module is required. ' +
                      'Use pip install pybloomfiltermmap3')

try:
    from progressbar import ProgressBar
    from progressbar import Bar
    from progressbar import Timer
    from progressbar import AdaptiveETA
    from progressbar import FormatCustomText
    from progressbar import Percentage
    from progressbar import SimpleProgress
except ImportError:
    raise ImportError('progressbar2 module is required. ' +
                      'Use pip install progressbar2')

from classes.kmerreader import KmerReader


class BFCounter():
    """
    Class for implementing Bloom Filter k-mer Counting algorithm

    Creates a Bloom Filter and check the k-mer is previously encountered or
    not. Only previously encountered k-mers are added to the Hash Table, which
    drastically reduced the size of the Hash Table.

    Referenced Paper:
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333
    """
    def __init__(self, reader, args=None):
        """
        :param  reader: KmerReader object
        :param  args: dictionary of additional arguments
        """
        if not isinstance(reader, KmerReader):
            raise TypeError('Reader should be of type KmerReader')
        # Default parameters
        self._error_rate = 1e-3
        self._verbose = False
        if args is None:
            pass
        elif isinstance(args, dict):
            if 'error_rate' in args:
                error_rate = args['error_rate']
                self._error_rate = error_rate
            if 'verbose' in args:
                verbose = args['verbose']
                self._verbose = verbose
        else:
            raise TypeError('Arguments should be of type dict')
        self._reader = reader
        self._heap = []  # min-heap for finding n most frequent items
        self._kmer_counter = dict()

    def nfrequent(self, n):
        """
        Returns a list of n most frequent k-mers and their counts

        :param  n: Desired number of most frequent k-mers to be returned
        :return:    List of n most frequent k-mers and their counts
        """
        self._count()
        self._heapify(n)
        return heapq.nlargest(n, self._heap)

    def _count(self):
        """
        Implements the Bloom Filter k-mer counting algorithm
        """
        # initialize Bloom Filter
        bf = BloomFilter(self._reader.total_kmer, self._error_rate, 'kmer_bf')
        if self._verbose:
            # initialize progress bar
            current = 0
            update_threshold = int(self._reader.total_kmer / 100)
            format_custom_text = FormatCustomText(
                'Hash Size: %(value).1f MB',
                dict(
                        value=0
                    )
            )
            print('Hashing...')
            bar = ProgressBar(max_value=self._reader.total_kmer, widgets=[
                Percentage(),
                ' ',
                SimpleProgress(format='(%s)' % SimpleProgress.DEFAULT_FORMAT,),
                ' ',
                Bar(),
                ' ',
                Timer(),
                ' ',
                AdaptiveETA(),
                ' ',
                format_custom_text
            ])
            bar.start()
        for kmer in self._reader.kmer():
            if kmer not in bf:  # not in Bloom Filter
                bf.add(kmer)
            else:  # in Bloom Filter
                if kmer in self._kmer_counter:  # in Hash Table
                    self._kmer_counter[kmer] += 1  # Increment
                else:
                    self._kmer_counter[kmer] = 2  # Add to Hash Table
            if self._verbose:
                # update progress bar
                current += 1
                if update_threshold == 0 or current % update_threshold == 0:
                    size = sys.getsizeof(self._kmer_counter) / (1024 ** 2)
                    bar.update(current,
                               format_custom_text.update_mapping(value=size))
        os.remove('kmer_bf')  # remove Bloom Filter from disk
        if self._verbose:
            bar.finish()
            print('Hashing Done!')

    def _heapify(self, n):
        # populate heap
        for i in range(n):
            self._heap.append((0, ''))
        for kmer, count in self._kmer_counter.items():
            if count > self._heap[0][0]:  # item is bigger than minimum item
                # replace minimum item with the recent one
                heapq.heappushpop(self._heap, (count, kmer))
