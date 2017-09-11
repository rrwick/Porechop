<p align="center"><img src="misc/porechop_logo_knife.png" alt="Porechop" width="600"></p>

You are looking at a strange, experimental branch of Porechop I made to split old 2D reads based on the hairpin adapter. If you're not here on purpose, I suggest you go back to the [master branch](https://github.com/rrwick/Porechop).

I hacked the code in this branch to do one thing only: split 2D reads. It skips all adapter searching and end trimming, jumping right to the middle adapter search using only the hairpin.

Usage:
`porechop -i input.fastq.gz -o output.fastq.gz --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0`

Like I said, it's experimental, so use at your own risk! :smile:
