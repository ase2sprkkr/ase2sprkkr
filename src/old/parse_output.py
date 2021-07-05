import numpy as np
import matplotlib.pyplot as plt

def get_scrape_info():
    info = [
        ('OUT_SCF', scrape_output),
    ]

    return info

def _skip_lines(fd, num):
    for ii in range(num):
        line = next(fd)
    return line

def scrape_output(filename, rdata=None):
    out = {
        'it' : [],
        'ERR' : [],
        'ETOT' : [],
        'EF' : [],
        'M' : [],
        'converged' : [],
        'E_band' : [],
    }
    ebs = {}
    with open(filename, 'r') as fd:
        for line in fd:
            if 'ERR' in line:
                items = line.split()
                out['it'].append(int(items[0]))
                out['ERR'].append(float(items[2]))
                out['EF'].append(float(items[5]))
                out['M'].append((float(items[10]), float(items[11])))
                items = _skip_lines(fd, 1).split()
                out['ETOT'].append(float(items[1]))
                flag = items[5] == 'converged'
                out['converged'].append(flag)

                out['E_band'].append(ebs)
                ebs = {}

                line = next(fd)
            elif 'SPRKKR-run for:' in line:
                run = line.replace('SPRKKR-run for:', '').strip()
                out['run'] = run

            elif ' E= ' in line:
                atom = line.split()[-1]
                line = _skip_lines(fd, 6)
                ebs[atom] = float(line.split()[1])

    return out

def get_plugin_info():
    info = [plot_err_etot, show_figures]

    return info

def plot_err_etot(df, data=None):
    fig, axs = plt.subplots(2, sharex=True)
    for ax in axs:
        ax.set_prop_cycle(
            plt.cycler(color=plt.cm.viridis(np.linspace(0, 1, len(df))))
        )

    vals = np.array([ii for ii in df['ERR']])
    ax = axs[0]
    for ir, val in enumerate(vals):
        print(ir, val)
        ax.semilogy(val, 'o', ls='None', mfc='None', mew=2, ms=6,
                    label=f'{ir}')
    ax.set_ylabel('ERR')

    vals = np.array([ii for ii in df['ETOT']])
    ax = axs[1]
    for ir, val in enumerate(vals):
        print(ir, val)
        ax.plot(val, 'o', ls='None', mfc='None', mew=2, ms=6,
                label=f'{ir}')
    ax.set_xlabel('IT')
    ax.set_ylabel('ETOT')
    fig.tight_layout()
    fig.savefig('err_etot.png')

def show_figures(df, data=None):
    plt.show()
    return data
