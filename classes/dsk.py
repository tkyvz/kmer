import math
import heapq
import os
import sys

try:
    from pybloomfilter import BloomFilter
except ImportError:
    raise ImportError('pybloomfilter module is required. ' +
                      'Use pip install pybloomfiltermmap3')

try:
    import mmh3
except ImportError:
    raise ImportError('mmh3 module is required. Use pip install mmh3')

try:
    from progressbar import ProgressBar
    from progressbar import UnknownLength
except ImportError:
    raise ImportError('progressbar2 module is required. ' +
                      'Use pip install progressbar2')

from classes.kmerreader import KmerReader


class DSK():
    """
    Class for implementing DSK, k-mer counting with very low memory algorithm

    Hashes each kmer and puts them into different files according to their
    hash values. By using target disk and memory spaces, determines how many
    distinct files should be used and in how many iterations the program needs
    to perform.

    Referenced Paper:
    http://minia.genouest.org/dsk/
    """
    def __init__(self, reader, args=None):
        """
        :param  reader: KmerReader object
        :param  args: dictionary of additional arguments
        """
        if not isinstance(reader, KmerReader):
            raise TypeError('Reader should be of type KmerReader')
        # Default parameters
        D = 50 * (1024 ** 3) * 8  # 50 GB in bits (target disk space)
        M = 4 * (1024 ** 3) * 8  # 4 GB in bits (target memory)
        self._error_rate = 1e-2  # Bloom Filter error rate
        self._verbose = False  # Verbose
        if args is None:
            pass
        elif isinstance(args, dict):
            if 'target_disk' in args:  # target disk space in GB
                _D = args['target_disk']
                D = _D * (1024 ** 3) * 8
            if 'target_memory' in args:  # target memory in GB
                _M = args['target_memory']
                M = _M * (1024 ** 3) * 8
            if 'error_rate' in args:  # error rate for Bloom Filter
                error_rate = args['error_rate']
                self._error_rate = error_rate
            if 'verbose' in args:
                verbose = args['verbose']
                self._verbose = verbose
        else:
            raise TypeError('Arguments should be of type dict')
        self._reader = reader
        b = (self._reader.k + sys.getsizeof('')) * 8  # size of kmer in bits
        b_disk = (self._reader.k + 1) * 8  # size of kmer and new line in file
        # calculate max unique kmers
        if (4 ** self._reader.k) > self._reader.total_kmer:
            v = reader.total_kmer
        else:  # cannot have more than 4^k kmer
            v = 4 ** self._reader.k
        self._niter = math.ceil(v * b_disk / D)
        self._np = math.ceil((v * b) / (0.7 * self._niter * M))
        self._capacity = v / (self._niter * self._np)
        self._heap = []
        self._efficient = 0.7 * M < b * v

    def should_use(self):
        """
        Determines whether DSK or BFCounter algorithm should be used based on
        the criterion provided in the referenced paper.
        :return: True if DSK should be used, False if BFCounter should be used
        """
        return self._efficient

    def nfrequent(self, n):
        """
        Returns a list of n most frequent k-mers and their counts

        :param  n: Desired number of most frequent k-mers to be returned
        :return:    List of n most frequent k-mers and their counts
        """
        if self._verbose:
            # Print number of iterations and partitions
            print('# of iterations:{}\t# of partitions: {}'.format(self._niter,
                                                                   self._np))
        self._populate(n)
        for i in range(self._niter):
            if self._verbose:
                print('Iteration #{}'.format(i + 1))
            self._write_files_for_iteration(i)  # perform iteration
            self._read_files_and_count()  # count files
            if self._verbose:
                print('Iteration #{} has been completed'.format(i + 1))
        return heapq.nlargest(n, self._heap)

    def _write_files_for_iteration(self, iter_no):
        """
        Performs one iteration of the DSK algorithm

        :param  iter_no: Index of the iteration to be performed
        """
        if self._verbose:
            # initialize progress bar
            print('Writing to files...')
            bar = ProgressBar(max_value=UnknownLength)
            bar.start()
            count = 0
        files = self._create_files()  # create files for this iteration
        for kmer in self._reader.kmer():
            h = mmh3.hash(kmer)  # hash the k-mer
            if h % self._niter == iter_no:  # belongs to this iteration
                j = int((h / self._niter) % self._np)  # determine partition
                files[j].write(kmer + '\n')  # write to file
                if self._verbose:
                    # update progress bar
                    count += 1
                    bar.update(count)
        for f in files:
            f.close()
        if self._verbose:
            bar.finish()
            print('Writing to files has been completed')

    def _create_files(self):
        """
        Creates files for each partition for a specific iteration
        """
        return [open('{}'.format(j), 'w') for j in range(self._np)]

    def _read_files_and_count(self):
        if self._verbose:
            print('Reading from files...')
        for j in range(self._np):
            if self._verbose:
                # initialize progress bar
                print('Partition #{}'.format(j + 1))
                bar = ProgressBar(max_value=UnknownLength)
                bar.start()
                count = 0
            bf = BloomFilter(
                self._capacity,
                self._error_rate,
                'kmer_bf'
            )
            kmer_counter = dict()
            with open(str(j), 'r') as f:  # open file for the current partition
                for kmer in f:
                    if kmer not in bf:  # not in Bloom Filter
                        bf.add(kmer)
                    else:  # in Bloom Filter
                        try:
                            kmer_counter[kmer] += 1  # in Hash Table
                        except KeyError:  # not in Hash Table
                            kmer_counter[kmer] = 2  # Add to Hash Table
                    if self._verbose:
                        # update progress bar
                        count += 1
                        bar.update(count)
            if self._verbose:
                bar.finish()
                print('Populating the heap...')
            for kmer, count in kmer_counter.items():
                if count > self._heap[0][0]:  # item is bigger than minimum
                    # replace minimum item with the recent one
                    # kmer.rstrip() is used to eliminate the new line
                    heapq.heappushpop(self._heap, (count, kmer.rstrip()))
            if self._verbose:
                print('Heap is populated')
                print(('Partition #{} has been completed with {:.1f} MB hash '
                       + 'table').format(
                       j + 1,
                       sys.getsizeof(kmer_counter) / (1024 ** 2)))
            os.remove(str(j))  # remove the partition file
            os.remove('kmer_bf')

    def _populate(self, n):
        """
        Populates min-heap
        :param  n: The number of elements the heap will have
        """
        for i in range(n):
            self._heap.append((0, ''))
