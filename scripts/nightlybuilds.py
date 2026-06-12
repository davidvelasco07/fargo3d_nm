import subprocess
import time

out=subprocess.check_output('ls -1 test_suite/*.py', shell=True).split()
TestList = []
for i in out:
    command = i.replace('test_suite/','test').replace('.py','')
    TestList.append (command)
    subprocess.call (["make", "mrproper"])
    subprocess.call (["make", command])

time.sleep (3)
out2=subprocess.check_output('grep "^test_suite/" tests.log', shell=True).split()
TestPerformed = []
for i in out2:
    if 'test_suite' in i:
        string = i.replace('test_suite/', 'test').replace('.py','')
        TestPerformed.append (string)


SetTest          = set(TestList)
SetTestPerformed = set(TestPerformed)

if len(SetTest-SetTestPerformed) > 0:
    with open('tests.log','a') as f:
        f.write ("FAILED. These tests have NOT been performed:\n")
        for test in SetTest-SetTestPerformed:
            f.write ("\t"+test)
            f.write ("\n")
