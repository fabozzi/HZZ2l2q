import os, commands
import string

comEnergy = '8TeV'

hmasses = ['400', '500', '600', '700', '800', '900', '1000' ]

templatecfg = "Higgs2l2qSkimPlusNtuple_Hpowheg_cfg_53.py" 

for masshyp in hmasses:
    file = open(templatecfg, "r")
    newfile = open('Higgs2l2qSkimPlusNtuple_Hpowheg_cfg_53'+'_'+masshyp+'.py', "w")
    print "new file: ", 'Higgs2l2qSkimPlusNtuple_Hpowheg_cfg_53'+'_'+masshyp+'.py'
    while 1:
        line = file.readline()
        newLine = line
#        if line == "":
        if len(line) == 0:
            break;
        if line.startswith("hMassHyp"):
            print line
            newLine = newLine.replace("700",masshyp)
            print newLine
        newfile.write(newLine)


    file.close()
    newfile.close()



