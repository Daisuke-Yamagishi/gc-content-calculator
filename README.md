# gc-content-calculator
This python script automatically detect your fasta files from your input directory and calculate gc content.

## How to use
You need to install the requirements.txt\
`pip install -r requirements.txt`\
Then, run the script.\
`python3 gc_content.py -i PATH_TO_YOUR_DIRECTORY_CONTAINING_FASTA_FILES' -o PATH_TO_OUTPUT_DIRECTORY`

## Output files
.tsv files and .png files (histogram)\
**Example (~.tsv)**
|ID|LENGTH|GC-count|GC_content|
|-|-|-|-|
|sequence1|100|12|1.2|
|sequence2|132|99|75|
|sequence3|678|87|8.895705521|\
