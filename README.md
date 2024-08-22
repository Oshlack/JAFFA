[![Docker Image Version](https://img.shields.io/docker/v/davidsongroup/jaffa?logo=docker&label=Docker%20Hub%20image&link=https%3A%2F%2Fhub.docker.com%2Fr%2Fdavidsongroup%2Fjaffa)](https://hub.docker.com/r/davidsongroup/jaffa)

JAFFA is a multi-step pipeline that takes either raw RNA-Seq reads, or pre-assembled transcripts, then searches for gene fusions. It will output the names and locations of candidate gene fusions along with the cDNA sequence of their breakpoints. JAFFA is based on the idea of comparing a transcriptome (e.g. in a cancer sample) against a reference transcriptome. In this way, it is a transcript-centric approach rather than a genome-centric approach like other fusion finders. In validation studies, JAFFA performed well over a range of read lengths - from 50bp to full-length transcripts and on single and paired-end reads.

See our [wiki](https://github.com/Oshlack/JAFFA/wiki) for downloads and instructions.

