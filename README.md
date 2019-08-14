# tmap
![tmap example layout](https://raw.githubusercontent.com/reymond-group/tmap/master/tmap/examples/mnist/mnist.jpg)

## Examples


## Getting started
tmap is a very fast visualization library for large, high-dimensional data sets. Currently, tmap is available for Python. tmaps graph layouts are based on the [OGDF](https://ogdf.uos.de/) library.

### Availability
| Language | Operating System | Status                 |
| -------- | ---------------- | ---------------------- |
| Python   | Linux            | Available              |
|          | Windows          | Available<sup>1</sup>  |
|          | macOS            | Unvailable<sup>2</sup> |
| R        |                  | Unvailable<sup>3</sup> |

<span class="small"><sup>1</sup>Works with
[WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)</span>  
<span class="small"><sup>2</sup>Availble by 19.08.2019</span>  
<span class="small"><sup>3</sup>FOSS R developers
[wanted](https://github.com/reymond-group/tmap)\!</span>

### Installation
tmap is installed using the conda package manager. Don't have conda? Download miniconda.

```bash
conda install -c tmap tmap
```

We suggest using faerun to plot the data layed out by tmap. But you can of course also use matplotlib (which might be to slow for large data sets and doesn't provide interactive features).

```bash
pip install faerun
# pip install matplotlib
```
### Tutorial and Documentation
See <a href="http://tmap.gdb.tools">http://tmap.gdb.tools</a>
