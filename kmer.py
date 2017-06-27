import time
import argparse


from classes.kmerreader import KmerReader
from classes.bfcounter import BFCounter
from classes.dsk import DSK


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


if __name__ == '__main__':
    # parser
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
                        default=1e-3,
                        type=check_between_zero_one)
    # DSK target disk space
    parser.add_argument('-d',
                        '--target-disk',
                        help='Target Disk Space in GB for DSK algorithm',
                        default=25,
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

    # Some algorithm parsing
    cli_args = parser.parse_args()
    verbose = cli_args.verbose
    quick = cli_args.quick
    args = dict()
    args['verbose'] = verbose
    args['error_rate'] = cli_args.error_rate
    args['target_disk'] = cli_args.target_disk
    args['target_memory'] = cli_args.target_memory
    n = cli_args.most_frequent
    k = cli_args.kmer_size
    file_name = cli_args.file_name

    # Initialize KmerReader object with the specified file and kmer size
    if verbose:
        start = time.time()

    reader = KmerReader(file_name, k)

    if verbose:
        end = time.time()
        print('{} {}-mers in {} counted in {:.2f} seconds'.format(
            reader.total_kmer, k, file_name, (end - start)
        ))

    if verbose:
        start = time.time()

    counter = DSK(reader, args)
    if not counter.should_use():  # if kmer_memory < 0.7 * target_memory
        counter = BFCounter(reader, args)
        if verbose:
            print('Selected algorithm: BFCounter')
    elif verbose:
        print('Selected algorithm: DSK')

    for item in counter.nfrequent(n):  # Get the n most frequent items
        count, kmer = item
        print('{}: {}'.format(kmer, count))

    if verbose:
        end = time.time()
        print('Total duration: {:.2f} seconds'.format((end - start)))
