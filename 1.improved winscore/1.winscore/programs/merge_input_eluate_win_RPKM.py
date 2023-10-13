import sys

inp=open(sys.argv[1])
ip=open(sys.argv[2])
output=open(sys.argv[3],'w')
output.write('\t'.join(["Transcript_id","GeneSymbol", "Input_rpkm", "IP_rpkm"])+'\n')

a=inp.readline()
b=ip.readline()

try:
    while True:
        a=inp.readline().rstrip('\n\r')
        b=ip.readline().rstrip('\n\r')
        if (a=='' and b==''):
            break
        if (a.split('\t')[0]!=b.split('\t')[0]):
            print ['Line match error: ', a, b]
        c=b.split('\t')[2]
        output.write(a+'\t'+c+'\n')
except EOFError as e:
    print e


    

