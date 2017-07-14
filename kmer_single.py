import argparse
from collections import defaultdict
import heapq
import math
import os
import sys
import time


try:
    from pybloomfilter import BloomFilter
except ImportError:
    raise ImportError('pybloomfilter module is required. ' +
                      'Use pip install pybloomfiltermmap3')

try:
    import mmh3
except ImportError:
    raise ImportError('mmh3 module is required. Use pip install mmh3')


def check_positive(value):
    """
    Checks value is a positive integer, returns True if so, else raise error.
    :param  value: value to be checked
    """
    try:
        ivalue = int(value)
        if ivalue <= 0:
            # is int but non-positive
            raise argparse.ArgumentTypeError(
                '{} is an invalid positive integer value'.format(value))
        return ivalue
    except ValueError:
        # not int
        raise argparse.ArgumentTypeError('{} is not an integer'.format(value))


def check_between_zero_one(value):
    """
    Checks value is between 0 and 1, returns True if so, else raise error.
    :param  value: value to be checked
    """
    try:
        fvalue = float(value)
        if 0 <= fvalue < 1:
            return fvalue
        # is float but not between 0 and 1
        raise argparse.ArgumentTypeError(
            'should be between 0 and 1'.format(value))
    except ValueError:
        # not float
        raise argparse.ArgumentTypeError(
            '{} is not a floating point number'.format(value))


def count_kmers(file_name, k, verbose=False):
    """
    Counts how many k-mers exists in a given file
    :param  file_name: Fastq file to be counted
    :param  k: K-mer size
    """
    if verbose:
        start = time.time()
        print('Counting kmers in {}'.format(file_name))
    total_kmers = 0
    with open(file_name, 'r') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # dna sequence
                total_kmers += len(line) - k  # eliminate new-line
            line_num += 1
    if verbose:
        end = time.time()
        print('{} kmers are counted in {:.2f} seconds'.format(
            total_kmers, end - start))
    return total_kmers


def parameters(total_kmers, target_disk, target_memory, k, verbose=False):
    """
    Calculates paramters
    :param  total_kmers: Number of total kmers in the file
    :param  target_disk: Target disk space
    :param  target_memory: Target memory
    :param  k: K-mer size
    :return number_of_iterations: Number of iterations for DSK
            number_of_partitions: Number of partitions for DSK
            bloom_filter_capacity: Bloom Filter Capacity
            dsk: True if DSK should be used, False if BFCounter should be used
    """
    if verbose:
        print('Calculating paramters...')
        print('Total K-mers: {}'.format(total_kmers))
        print('Target Disk: {:.2f}GB'.format(target_disk / 1024**3 / 8))
        print('Target Memory: {:.2f}GB'.format(target_memory / 1024**3 / 8))

    b = (k + sys.getsizeof('')) * 8  # size of kmer in bits
    b_disk = (k + 1) * 8  # size of kmer in disk in bits
    # max unique kmers
    if (4 ** k) < total_kmers:
        v = 4 ** k
    else:
        v = total_kmers
    # DSK Paramters
    number_of_iterations = int(math.ceil(total_kmers * b_disk / target_disk))
    number_of_partitions = int(
        math.ceil((v * b) / (0.7 * target_memory * number_of_iterations)))
    use_dsk = 0.7 * target_memory < b * v
    # Bloom Filter Capacity
    if dsk:
        bf_capacity = total_kmers / (number_of_iterations *
                                     number_of_partitions)
    else:
        bf_capacity = total_kmers

    if verbose:
        print('Parameters are calculated')
        if use_dsk:
            print('Algorithm: DSK')
            print('# of iterations: {}'.format(number_of_iterations))
            print('# of partitions: {}'.format(number_of_partitions))
        else:
            print('Algorithm: BFCounter')
        print('Bloom Filter Capacity: {}'.format(bf_capacity))

    return number_of_iterations, number_of_partitions, bf_capacity, use_dsk


def dsk(file_name, k, n, capacity, error_rate, iters, parts, verbose=False):
    """
    Implementation of DSK, k-mer counting with very low memory algorithm

    Hashes each kmer and puts them into different files according to their
    hash values. By using target disk and memory spaces, determines how many
    distinct files should be used and in how many iterations the program needs
    to perform.

    Referenced Paper:
    http://minia.genouest.org/dsk/
    """
    if verbose:
        start = time.time()

    # Assign functions to local variables for performance improvement
    hash_function = mmh3.hash
    heap_pushpop = heapq.heappushpop

    CHUNK_LIMIT = math.floor(capacity / 10)  # write approximately in 10 calls

    heap = []
    for i in range(n):
        heap.append((0, ''))
    for it in range(iters):  # iteration
        if verbose:
            start_iter = time.time()
            print('Iteration#{} started.'.format(it + 1))
        files = [open('{}'.format(j), 'w') for j in range(parts)]  # open files

        # Write to files in chunks to have less file.write calls
        chunks = [[] for j in range(parts)]

        # Assign functions to local variables for performance improvement
        writers = [files[j].write for j in range(parts)]
        chunk_appender = [chunks[j].append for j in range(parts)]
        chunk_cleaner = [chunks[j].clear for j in range(parts)]
        chunk_joiner = ''.join

        with open(file_name, 'r') as f:
            line_num = 0
            for line in f:
                if line_num % 4 == 1:  # dna sequence
                    kmer_count = len(line) - k
                    for i in range(kmer_count):
                        kmer = line[i:i + k]
                        h = hash_function(kmer)
                        if h % iters == it:  # belongs to this iteration
                            j = (h / iters) % parts
                            _j = int(j)
                            chunk_appender[_j](kmer + '\n')
                            if len(chunks[_j]) == CHUNK_LIMIT:
                                # write to file
                                writers[_j](chunk_joiner(chunks[_j]))
                                chunk_cleaner[_j]()
                line_num += 1

        # Write remaining kmers
        for j in range(parts):
            writers[j](chunk_joiner(chunks[j]))

        for f in files:
            f.close()  # close files

        del chunks

        if verbose:
            end_disk_write = time.time()
            print('Disk write is completed in {:.2f} seconds.'.format(
                end_disk_write - start_iter
            ))

        for j in range(parts):
            bf = BloomFilter(capacity, error_rate, 'kmer_bf')

            kmer_counter = defaultdict(lambda: 1)

            # Assign functions to local variables for performance improvement
            add_to_bf = bf.add

            if verbose:
                start_partition = time.time()
                print('Reading partition#{} started.'.format(j + 1))

            with open(str(j), 'r') as f:
                for kmer in f:
                    if kmer not in bf:  # not in Bloom Filter
                        add_to_bf(kmer)
                    else:  # in Bloom Filter
                        kmer_counter[kmer] += 1

            if verbose:
                end_partition = time.time()
                print('Reading partition#{} is completed '.format(j + 1) +
                      'in {:.2f} seconds.'.format(
                      end_partition - start_partition))
                start_populate = time.time()
                print('Populating the heap...')

            for kmer, count in kmer_counter.items():
                # insert to the heap if count is bigger than minimum
                if count > heap[0][0]:
                    heap_pushpop(heap, (count, kmer.rstrip()))

            if verbose:
                end_populate = time.time()
                print('Heap is populated in {:.2f} seconds.'.format(
                    end_populate - start_populate
                ))

            os.remove(str(j))
            os.remove('kmer_bf')

    if verbose:
        end = time.time()
        print('DSK Duration: {:.2f} seconds.'.format(end - start))
    return heap


def bf_counter(file_name, k, n, capacity, error_rate, verbose=False):
    """
    Implementation of Bloom Filter k-mer Counting algorithm

    Creates a Bloom Filter and check the k-mer is previously encountered or
    not. Only previously encountered k-mers are added to the Hash Table, which
    drastically reduced the size of the Hash Table.

    Referenced Paper:
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333
    """
    if verbose:
        start = time.time()
        print('BFCounter started.')

    heap = []
    for i in range(n):
        heap.append((0, ''))

    bf = BloomFilter(capacity, error_rate, 'kmer_bf')

    kmer_counter = defaultdict(lambda: 1)

    # Assign functions to local variables for performance improvement
    add_to_bf = bf.add
    heap_pushpop = heapq.heappushpop

    with open(file_name, 'r') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # dna sequence
                kmer_count = len(line) - k
                for i in range(kmer_count):
                    kmer = line[i:i + k]
                    if kmer not in bf:  # not in Bloom Filter
                        add_to_bf(kmer)
                    else:  # in Bloom Filter
                        kmer_counter[kmer] += 1
            line_num += 1
    if verbose:
        end_hash = time.time()
        hash_table_size = sys.getsizeof(kmer_counter) / (1024 ** 2)
        print('Hash table is created in {:.2f} seconds.'.format(
            end_hash - start))
        print('Hash table size: {:.2f} MB.'.format(hash_table_size))
        start_populate = time.time()
        print('Populating the heap...')

    for count, kmer in kmer_counter.items():
        # insert to the heap if count is bigger than minimum
        if count > heap[0][0]:
            heap_pushpop(heap, (count, kmer))

    if verbose:
        end_populate = time.time()
        print('Heap is populated in {:.2f} seconds.'.format(
            end_populate - start_populate
        ))

    os.remove('kmer_bf')
    if verbose:
        end = time.time()
        print('BFCounter is completed in {:.2f} seconds.'.format(end - start))

    return heap


if __name__ == '__main__':
    desc = 'A program for counting most frequent k-mers in .fastq files.'
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.RawTextHelpFormatter)
    # file name
    parser.add_argument('-f',
                        '--file-name',
                        required=True,
                        help='Name of the .fastq file to be processed',
                        type=str)
    # kmer size
    parser.add_argument('-k',
                        '--kmer-size',
                        required=True,
                        help='Length of k-mers',
                        type=check_positive)
    # most frequent count
    parser.add_argument('-n',
                        '--most-frequent',
                        required=True,
                        help='Number of most frequent k-mers to be outputted',
                        type=check_positive)
    # Bloom Filter error rate for BFCounter
    parser.add_argument('-e',
                        '--error-rate',
                        help='Bloom Filter error rate for BFCounter algorithm',
                        default=1e-2,
                        type=check_between_zero_one)
    # DSK target disk space
    parser.add_argument('-d',
                        '--target-disk',
                        help='Target Disk Space in GB for DSK algorithm',
                        default=50,
                        type=check_positive)
    # DSK target disk memory
    parser.add_argument('-m',
                        '--target-memory',
                        help='Target Memory in GB for DSK algorithm',
                        default=4,
                        type=check_positive)
    # Verbose
    parser.add_argument('-v',
                        '--verbose',
                        help='Verbose',
                        action='store_true')

    # Parameters
    cli_args = parser.parse_args()
    verbose = cli_args.verbose  # verbose
    error_rate = cli_args.error_rate  # bloom filter error rate
    target_disk = cli_args.target_disk * (1024**3) * 8  # target disk
    target_memory = cli_args.target_memory * (1024**3) * 8  # target memory
    n = cli_args.most_frequent  # n, the number of kmers to be returned
    k = cli_args.kmer_size  # k, size of kmer
    file_name = cli_args.file_name  # fastq file name

    start = time.time()

    # Count total k-mers
    total_kmers = count_kmers(file_name, k, verbose=verbose)
    # Calculate paramters
    iters, parts, capacity, is_dsk = parameters(total_kmers,
                                                target_disk,
                                                target_memory,
                                                k,
                                                verbose=verbose)
    if is_dsk:  # DSK Algorithm implementation
        heap = dsk(
            file_name, k, n, capacity, error_rate, iters, parts, verbose
        )
    else:  # BFCounter Algorithm implementation
        heap = bf_counter(
            file_name, k, n, capacity, error_rate, verbose
        )

    for count, kmer in heapq.nlargest(n, heap):
        print('{}: {}'.format(kmer, count))

    end = time.time()
    print('Duration: {:.2f} seconds'.format(end - start))
