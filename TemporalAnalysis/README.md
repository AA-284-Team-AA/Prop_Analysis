# Temporal performance analysis
## Usage

To use the code, use
```bash
python time_analysis.py <input_file_name>
```
By default, the code will look for `input.ini`. Instructions on what to provide to the input file are written as comments and should be self-explanatory.

To modify the `mdot_ox` profile, change the implementation of `get_desired_mdot_ox` in `time_analysis.py`.

The code uses python 3. I don't remember what packages I installed myself. If you want to work in the same environment as me, start with a fresh python 3 install (using `virtualenv` for instance) and do
```bash
pip install -r my_env.txt
```