META            = 'Hs_meta.csv'

## non default parameters
ANNOTATION      = '/annotation/Homo_sapiens.GRCh38.97.gtf' 
GENOME_FASTA    = '/indexes/bowtie/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
GENEPRED        = '/annotation/Homo_sapiens.GRCh38.97.genePred.wgn'

GENOME_INDEX    = "/indexes/hisat2/Homo_sapiens.GRCh38.dna.primary_assembly"
SEGEMEHL_INDEX  = "/indexes/segemehl/Homo_sapiens.GRCh38.dna.primary_assembly.idx"
BWA_INDEX       = "/indexes/bwa/Homo_sapiens.GRCh38.dna.primary_assembly"
BOWTIE2_INDEX   = "/indexes/bowtie2/Homo_sapiens.GRCh38.dna.primary_assembly"
STAR_INDEX      = "/indexes/star/Homo_sapiens.GRCh38.dna.primary_assembly"
BOWTIE_INDEX    = "/indexes/bowtie/Homo_sapiens.GRCh38.dna.primary_assembly"

CPUS = "12"

PREPROCESSOR    = "trimmomatic"
PREPROCESSOR_PARAMS = "MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:26 CROP:150"

FIX_READ_HEADER = 'True'
SAM_SORT_MM = '3G'

