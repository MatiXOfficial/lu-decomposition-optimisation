import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

Y_MAX = 7000

RANGE = range(1, 10)

FILES = {
    "-": [f'results/lu{i}_results.txt' for i in RANGE],
    "-O2": [f'results/lu{i}_results_O2.txt' for i in RANGE],
}

NAMES = {
    "-": [f'lu{i}' for i in RANGE],
    "-O2": [f'lu{i} -O2' for i in RANGE],
}

all_datasets = {}
for option, files in FILES.items():
    datasets = []
    for file in files:
        dataset = pd.read_csv(file)
        dataset.index = dataset['Size']
        dataset.drop(columns=['Size'], inplace=True)
        
        sizes = np.array(dataset.index)
        n_Mflop = (sizes * (4 * sizes * sizes - 3 * sizes - 1)) / 6_000_000
        dataset = 1. / dataset
        dataset = dataset.multiply(n_Mflop, axis=0)
        dataset.columns = [f'MFLOPS_{i}' for i in range(len(dataset.columns))]

        dataset['mean'] = dataset.mean(axis=1)
        dataset['std'] = dataset.std(axis=1)
        datasets += [dataset]

    all_datasets[option] = datasets


for i, (option, datasets) in enumerate(all_datasets.items()):
    # one versus the other
    for j, dataset in enumerate(datasets):
        name = NAMES[option][j]
        x = np.array(dataset.index)
        mean = np.array(dataset['mean'])
        std = np.array(dataset['std'])

        if j == 0:
            name_0 = name
            x_0 = x
            mean_0 = mean
            std_0 = std
        else:
            name = NAMES[option][j]
            x = np.array(dataset.index)
            mean = np.array(dataset['mean'])
            std = np.array(dataset['std'])

            for cmp_name, cmp_x, cmp_mean, cmp_std in [(name_0, x_0, mean_0, std_0), (name_prev, x_prev, mean_prev, std_prev)]:
                plt.errorbar(x, mean, std, marker='.', linewidth=0.9, capsize=4, capthick=1, elinewidth=0.8, label=name)
                plt.errorbar(cmp_x, cmp_mean, cmp_std, marker='.', linewidth=0.9, capsize=4, capthick=1, elinewidth=0.8, label=cmp_name)
                plt.title(f'{name_prev} vs {name}')
                plt.xlabel('rozmiar')
                plt.ylabel('MFLOPS')
                plt.ylim(bottom=0, top=Y_MAX)
                plt.legend()
                plt.tight_layout()
                plt.savefig(f'plots/{option}/{cmp_name} vs {name}.jpg')
                # plt.show()
                plt.clf()

        name_prev = name
        x_prev = x
        mean_prev = mean
        std_prev = std

    # all in one
    for j, dataset in enumerate(datasets):
        name = NAMES[option][j]
        x = np.array(dataset.index)
        mean = np.array(dataset['mean'])
        std = np.array(dataset['std'])

        plt.plot(x, mean, marker='.', linewidth=0.9, label=name)

    plt.title(f'Porównanie optymalizacji ({option})')
    plt.xlabel('rozmiar')
    plt.ylabel('MFLOPS')
    plt.ylim(bottom=0, top=Y_MAX)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/{option}/all.jpg')
    # plt.show()
    plt.clf()
