import os

path = os.curdir
files = os.listdir(path)
fastq_Files = [i for i in files if i.endswith('.fastq')]


for fastq_file in fastq_Files:
    filename = fastq_file.split('.')[0]
    filename += '_noPlus.fastq'
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        print("File_before:" + fastq_file + " Length:" + str(len(lines)))
        for i in range(2, len(lines), 4):
            lines[i] = '+\n'
        print("File_after:" + fastq_file + " Length:" + str(len(lines)))
    with open(filename, 'w+') as out_file:
        out_file.writelines(i for i in lines)
