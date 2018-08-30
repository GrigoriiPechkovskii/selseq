'''intended for additional code'''
print('start selseq_extra')
from selseq_constant import *
from selseq_main import SequenceFasta

def nothing(obj):    
    return obj

def generator_fasta2(directory,type_fasta='faa',type_output=None,fun=None):
    '''
    '''
    #global fasta  
    files = os.listdir(directory)
    for file in files:
        if file.endswith(type_fasta):
            fasta = SequenceFasta(directory+file)
            fasta.seq_process(strip=True)
    for ind in range(fasta.seq_len):
        #yield ind
        yield fasta

def generator_fasta(directory,type_fasta='faa',fun=nothing,type_output=None):
    '''
    '''
    #global fasta  
    files = os.listdir(directory)
    for file in files:
        if file.endswith(type_fasta):
            fasta = SequenceFasta(directory+file)
            fasta.seq_process(strip=True)
            yield fun(fasta)

    

class generator_fasta3():

    def __init__(self,directory,type_fasta='faa'):
        self.directory = directory
        self.type_fasta = type_fasta

    def generator_index(self):        
        files = os.listdir(self.directory)
        for file in files:
            if file.endswith(self.type_fasta):
                fasta = SequenceFasta(self.directory+file)
                fasta.seq_process(strip=True)
                for ind in range(fasta.seq_len):
                    yield ind
                    #yield fasta

        

if __name__ == '__main__':
    for f in generator_fasta(REDATA_DIRECTORY):    
        print(f)

print('end selseq_extra')
