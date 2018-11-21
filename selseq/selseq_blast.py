'''for blast'''
from selseq_constant import *
import os.path
from selseq_main import blast
import os
print('start selseq_blast')

def db_for_blast(assemble_files,name):
    '''Make faa file from assemble for total blast'''
    
    with open(REDATA_DIRECTORY + name + '.faa','a') as join_assemble:
        for file in assemble_files:
            
            with open(REDATA_DIRECTORY+file,'r') as assemble_file_opened:
                assemble_file_read = assemble_file_opened.read()
                join_assemble.write(assemble_file_read)

    path_join_assemble = REDATA_DIRECTORY + name + '.faa'
    return path_join_assemble 
                #join_assemble.remove()

def blast_total():
    n = 0
    for index in range(len(ASSEMBLE_FILES)):    
        n+=1
        print(ASSEMBLE_FILES[0])    
        path_assemble = db_for_blast(ASSEMBLE_FILES,'join_assemble' + str(n))
        name_assemble = os.path.basename(path_assemble)
        #blast()
        blast(name_assemble,ASSEMBLE_FILES[0],name_assemble[0:-4] + '_db',name_assemble[0:-4] + '_tbl.csv')
        os.remove(REDATA_DIRECTORY + name_assemble[0:-4] + '_db.phr')
        os.remove(REDATA_DIRECTORY + name_assemble[0:-4] + '_db.pin')
        os.remove(REDATA_DIRECTORY +name_assemble[0:-4] + '_db.psq')
        os.remove(REDATA_DIRECTORY +name_assemble)
        ASSEMBLE_FILES.pop(0)

if __name__ == '__main__':
    blast_total()

print('end selseq_blast')
