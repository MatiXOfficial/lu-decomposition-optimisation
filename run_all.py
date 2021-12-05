import subprocess

FILES = [f'lu/lu{i}.c' for i in range(1, 10)]
OPTIONS = [("", ""), ("-O2", "O2")]

for option, option_string in OPTIONS:
    for file in FILES:
        print('-------------------------------')
        print(file + ' ' + option)
        subprocess.run([f"gcc {file} -mavx {option}"], shell=True)
        subprocess.run([f"./a.out {option_string}"], shell=True)
