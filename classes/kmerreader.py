import os.path


class KmerReader():
    """
    A file reader class for reading K-mers from FASTQ files.
    """
    def __init__(self, file_name, kmer_size):
        """
        :param  file_name: FASTQ file to be processed
        :param  kmer_size: kmer size to be counted in the file
        """
        if (os.path.isfile(file_name)):
            self._file = file_name
            self.k = kmer_size
            self.total_kmer = self._count()
            if (self.total_kmer == 0):  # invalid file
                raise TypeError(self._file + ' is not a valid FASTQ file.')
        else:
            raise FileNotFoundError('File ' + file_name + ' not found!')

    def kmer(self):
        """
        Iterates over kmers in a file
        """
        with open(self._file, 'r') as f:
            line_num = 0
            for line in f:
                if (line_num % 4 == 1):
                    kmers_in_line = len(line) - self.k  # eliminate new line
                    for i in range(kmers_in_line):
                        kmer = line[i:self.k + i]  # yield kmer for use
                        yield kmer
                line_num += 1
        return None

    def _count(self):
        """
        Counts the total number of kmers in the file.
        If the file is not in FASTQ format, returns 0.
        """
        count = 0
        with open(self._file, 'r') as f:
            line_num = 0
            for line in f:
                if (line_num % 4 == 0):  # sample id
                    if (not line.startswith('@')):
                        count = 0
                        break
                elif (line_num % 4 == 1):  # dna sequence
                    count += len(line) - self.k  # eliminate new line char
                elif (line_num % 4 == 2):  # comment
                    if (not line.startswith('+')):
                        count = 0
                        break
                else:  # quality (do nothing)
                    pass
                line_num += 1
        return count
