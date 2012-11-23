import os, commands
import string

comEnergy = '8TeV'

hmasses = ['400', '425', '450', '475', '500', '525', '550', '575', '600', '650', '700', '750', '800', '850', '900', '950', '1000']

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
            newLine = newLine.replace("400",masshyp)
            print newLine
        newfile.write(newLine)


    file.close()
    newfile.close()



