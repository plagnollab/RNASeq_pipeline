# Predict functional consequences of splice variant using SGSeq
library(data.table)
library(dplyr)
library(stringr)
library(SGSeq)
library(intervals)

# plot variants from sgv object
plotSGV <- function(mysgv){
  # bind together all elements of list as one object
  all <- unlist(mysgv)
  # separate by type
  exons <- SGSeq::type(all) == "E"
  colours <- ifelse( exons, "red", "blue")
  # create interval matrices
  mat <- as.matrix(ranges(all) )
  mat[,2] <- mat[,1] + mat[,2]
  mat <- Intervals(mat)
  transcript_names <- paste(SGSeq::txName(all), collapse = "\n" )
  transcript_names <- gsub("$^", "novel", transcript_names)
  row.names(mat) <-  transcript_names
  plot( mat, col = colours, use_names = TRUE, names_cex = 0.75, use_points = TRUE )
}

# testing
species <- "mouse"
sgseq_res <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Control_KO/Nicol_FUS_KO_Control_KO_res_clean_novel.tab"
support_frame <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/nicol_KO_SGSeq_support.tab"
code <- "Nicol_FUS_KO"
condition <- "condition_HOM"
sgf_object <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Nicol_FUS_KO_sgv_novel.RData"



option_list <- list(
  make_option(c('--support.tab'), help='', default = "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i_new_support.tab"),
  make_option(c('--code'), help='', default = "F210I_embryonic_brain_norm"),
  make_option(c('--case.condition'), help='', default = "HOM"), # now deprecated
  make_option(c('--sgf_object')),
  make_option(c('--sgseq_res'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"),
  make_option(c('--output.dir'), help='', default="")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sgf_object <- opt$sgf_object
support.tab <- opt$support.tab
code <- opt$code
gtf <- opt$gtf
sgseq_res <- opt$sgseq_res
outFolder <- opt$output.dir 
species <- opt$species 


# load in SGSeq R object
message("loading in SGSeq object - may take a while")
load(sgf_object)
# read in SGSeq results
d <- as.data.frame(fread(sgseq_res), stringsAsFactors=FALSE)

cassette_exons <- subset(d, variantType == "SE:I" & padj < 0.05)

# predict effect of cassette splicing variants
if( species == "mouse"){
  genome <- "BSgenome.Mmusculus.UCSC.mm10"
  transcripts <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
}
#TODO human, fly, worm

# load in genome and transcript list
library(genome, character.only = TRUE)
library(transcripts, character.only = TRUE)

txdb <- eval(parse(text = transcripts))
genome <- eval(parse(text = genome))


# subset sgv_novel for just the significant cassette exons (included)
cassette_ids <- unique(cassette_exons$groupID)
cassettes <- sgv_novel[ which(SGSeq::eventID(sgv_novel) %in% cassette_ids ) ]

message(paste("Number of cassette exons:", length(cassettes) ))

# for each cassette exon
# put both inclusion and exclusion event through predictVariantEffects
# if skiptic or cryptic, the novel variant won't be predicted
# but the annotated event will and will predict the effect of inclusion/exclusion
cassette_predictions <- list()
for( i in 1:length(cassette_ids) ){
  variants <- cassettes[ eventID(cassettes) == cassette_ids[i] ] 
  plotSGV(variants)
  vep <- predictVariantEffects(variants, tx = txdb, genome = genome, output = "full", cores = 2)
  # go through each row of prediction
  pred <- list()
  for( j in 1:nrow(vep)){
    ref <- vep$protein_ref_seq[j]
    var <- vep$protein_var_seq[j]
    # check if protein sequences are identical 
    ## if so the splice variant doesn't change protein as not in CDS
    ### so check cDNA sequence instead
    protein_change <- vep[j,]$protein_variant_type
    RNA_change <- vep[j,]$RNA_variant_type
    
    if( ref == var & vep[j,]$protein_variant_type == "no_change" ){
      protein_change <- "CDS"
      ref <- vep$RNA_ref_seq[j]
      var <- vep$RNA_var_seq[j]
    }
    # align shorter sequence against longer one to get sequence of inserted string
    if( str_length(var) > str_length(ref) ){
      p <- pairwiseAlignment(ref, var, type = "global")
    }else{
      p <- pairwiseAlignment(var, ref, type = "global")
    }
    
    segment <- unlist( deletion(p) )
    
    # if a frameshift occurs and a totally new protein seq is created then perform a local alignment and get out the novel protein seq
    if( length(segment@start) > 1 ){
      if( str_length(var) > str_length(ref) ){
        p <- pairwiseAlignment(ref, var, type = "local")
      }else{
        p <- pairwiseAlignment(var, ref, type = "local")
      }
      
      # get the protein tag at the end
      aligned <- as.character(pattern(p))
      segment_seq <- str_sub( var,  start = str_length(aligned) +1 )
    }else{
      # if a straightforward insertion then get out the central sequence
      segment_seq <- str_sub( ref, start = segment@start, end = segment@start + segment@width)
    }
    
    
    print( list( RNA_change = RNA_change, protein_change = protein_change, segment_seq = segment_seq))
    pred[[j]] <- list( RNA_change = RNA_change, protein_change = protein_change, segment_seq = segment_seq)
    
  }
  cassette_predictions[[i]] <- pred
}

# predict variant effects - a little slow
message("predicting splice variant outcomes on protein sequence")
vep <- predictVariantEffects(cassettes[1:6], tx = txdb, genome = genome, output = "full", cores = 2)

# vep is a dataframe - can be manipulated
vep <- predictVariantEffects(cassette_exon_variants[11:12], tx = txdb, genome = genome, output = "full", cores = 2)

# can use pairwiseAlignment function from biostrings to extract the protein sequence of the central exon
ref <- vep$protein_ref_seq[2]
var <- vep$protein_var_seq[2]
p <- pairwiseAlignment(var, ref, type = "global")
central_seq <- unlist(deletion(p))
# ref is the inclusion protein
central <- str_sub( ref, start = central_seq@start, end = central_seq@start + central_seq@width)








plotSGV(cassette_exon_variants[11:12])


MFSVRIVTADYYMASPLPGLDTCQSPLTQLPVKKVPVVRVFGATPAGQKTCLHLHGIFPYLYVPYDGYGQQPESYLSQMAFSIDRALNVALGNPSSTAQHVFKVSLVSGMPFYGYHEKERHFMKIYLYNPAMVKRICELLQSGAIMNKCYQPHEAHIPYLLQLFIDYNLYGMNLINLAAVKFRKARRKGNASHATGLFKHQLSGNSPAGTLFRWEEDEIPSSLLLEGVEPLSTCELEVDAVAADILNRLDIEAQIGGNPGLQAIWEDEKQRRRNRNESSQISQPESQDCRFVPATESEKQFQKRLQEVLKQNDFSVTLSGSVDYSNGSQEFSAELTLHSEILSPEMLPCSPANMIEVHKDTDLSKGNTKHKVEEALINEEAILNLIENSQTFQPLTQRLSETPVFMGSSPDESLVHLLAGLESDGYQGEKNRMPLPCHSFGESQNPQNSDDEENEPQIEKEEMELSVVMSQRWDSDIEEHCAKKRSLCRNAHRSSTEEDDSSSEEEMEWTDNSLLFANLSIPQLDGTADENSDNPLNNENSRAHSSVIATSKLSVRPSIFHKDAATLEPPSSAKITFQCKHTSALSSHVLNKDGLTEDLSQPNSTEKGRDNSVTFTKESTYSMKYSGSLSSTVHSDNSHKEICKKDKSLPVSSCESSVFDYEEDIPSVTRQVPSRKYSNMRKIEKDASCIHVNRHISETILGKNSFNFADLNHSKRKLSSEGNEKGNSTSLSGVFPSSLTENCDLLPSSGENRSMAHSLESITDESGLNKLKIRYEEFQEHKMEKPSLSQQAAHYMFFPSVVLSNCLTRPQKLSPVTYKLQSGNKPSRLKLNKKKLIGLQETSTKSTETGATKDSCTHNDLYTGASEKENGLSSDSAKATHGTFENKPPTEHFIDCHFGDGSLEAEQSFGLYGNKYTLRAKRKVNYETEDSESSFVTQNSKISLPHPMEIGENLDGTLKSRKRRKMSKKLPPVIIKYIIINRFRGRKNMLVKLGKIDSKEKQVILTEEKMELYKKLAPLKDFWPKVPDSPATKYPIYPLTPKKSHRRKSKHKSAKKKPGKQHRTNSENIKRTLSFRKKRTHAVLSPPSPSYIAETEDCDLSYSDVMSKLGFLSERSTSPINSSPPRCWSPTDPRAEEIMAAAEKESMLFKGPNVYNTKTVSPRVGKASRARAQVKKSKARLANSSVVTNKRNKRNQTTKLVDDGKKKPRAKQKQRANEKSLSRKHAIPADEKMKPHSEAELTPNHQSVSELTSSSGAQALSKQKEMSQTGPAVDHPLPPAQPTGISAQQRLSNCFSSFLESKKSVDLRTFPSSRDDSHSSVVYSSIGPGISKINIQRSHNQSAMFTRKETTLIQKSIFDLSNHLSQVAQSTQVCSGIISPKTEESSSTQKNCGSSMGKLNEYRSSLESKPEQVCAPNFLHCKDSQQQTVSVSEQSKTSETCSPGNAASEESQTPNCFVTSLKSPIKQIAWEQKQRGFILDMSNFKPEKVKQRSLSEAISQTKALSQCKNQNVSTPSVFGEGQSGLAVLKELLQKRQQKAQSTNVVQDSTSTHQPDKNISVSNEHKKANKRTRPVTSPRKPRTPRRTKPKEQTPRRLKVDPLNLQTSGHLDNSLSDDSPILFSDPGFESCYSLEDSLSPEHNYNFDINTIGQTGFCSFYSGSQFVPADQNLPQKFLSDAVQDLFPGQAIDKSELLSHDRQSCSEEKHHVSDSSPWIRASTLNPELFEKVAMDNNENHRHSQWKNSFHPLTSHSNSIMESFCVQQAENCLTEKSRLNRSSVSKEVFLSLPQANSSDWIQGHNRKEADQSLHSANTSFTTILSSPDGELVDAASEDLELYVSRNNDVLTPTPDSSPRSTSSPLQSKNGSFTPRTAHILKPLMSPPSREEIVATLLDHDLSEAIYQEPFCSNPSDVPEKPREIGGRLLMVETRLPNDLIEFEGDFSLEGLRLWKTAFSAMTQNPRPGSPLRNGQAVVNKESSNSHKMVEDKKIVIMPCKYAPSRQLVQAWLQAKEEYERSKKLPKTELTPVTKSAENVSPSLNPGDTCAVSPQVDKCPHTLSSSAHTKEEVSKSQIALQTSTTGCSQTLLAAASAAVPEEDEDDNDNCYVSYSSPDSPGIPPWQQAASPDFRSLNGDDRHSSPGKELCSLAVENFLKPIKDGIQKSSCSESWEPQVISPIHARARTGKWDPLCLHSTPVMQRKFLEKLPEATGLSPLSVEPKTQKLYNKKGSDADGLRRVLLTTQVENQFAAVNTPKKETSQIDGPSLNNTYGFKVSIQNLQEAKALHEIQNLTLISVELHARTRRDLQPDPEFDPICALFYCISSDTPLPDTEKTELTGVIVIDKDKTVTHQDIRSQTPLLIRSGITGLEVTYAADEKALFQEITNIIKRYDPDILLGYEIQMHSWGYLLQRAAALSVDLCQMISRVPDDKIENRFAAERDDYGSDTMSEINIVGRITLNLWRIMRNEVALTNYTFENVSFHVLHQRFPLFTFRVLSDWFDNKTDLYRWKMVDHYVSRVRGNLQMLEQLDLIGKTSEMARLFGIQFLHVLTRGSQYRVESMMLRIAKPMNYIPVTPSIQQRSQMRAPQCVPLIMEPESRFYSNSVLVLDFQSLYPSIVIAYNYCFSTCLGHVENLGKYDEFKFGCTSLRVPPDLLYQIRHDVTVSPNGVAFVKPSVRKGVLPRMLEEILKTRLMVKQSMKSYKQDRALSRMLNARQLGLKLIANVTFGYTAANFSGRMPCIEVGDSIVHKARETLERAIKLVNDTKKWGARVVYGDTDSMFVLLKGATKEQSFKIGQEIAEAVTATNPRPVKLKFEKVYLPCVLQTKKRYVGYMYETLDQKEPVFDAKGIETVRRDSCPAVSKILERSLKLLFETRDISLIKQYVQRQCMKLVEGKASIQDFIFAKEYRGSFSYRPGACVPALELTRKMLAYDRRSEPRVGERVPYVIIYGTPGLPLIQLIRRPAEVLQDPTLRLNATYYITKQILPPLARIFSLIGIDVFSWYQELPRIQKATSSSRSELEGRKGTISQYFTTLHCPVCDDLTQHGICSKCRSQPQHVAIILNQEIRELERKQEQLIKICRNCTGSFDRHIPCVSLNCPVLFKLSRVNRELSKAPYLRQLLDQF
MFSVRIVTADYYMASPLPGLDTCQSPLTQLPVKKVPVVRVFGATPAGQKTCLHLHGIFPYLYVPYDGYGQQPESYLSQMAFSIDRALNVALGNPSSTAQHVFKVSLVSGMPFYGYHEKERHFMKIYLYNPAMVKRICELLQSGAIMNKCYQPHEAHIPYLLQLFIDYNLYGMNLINLAAVKFRKARRKGNASHATGLFKHQLSGNSPAGTLFRWEEDEIPSSLLLEGVEPLSTCELEVDAVAADILNRLDIEAQIGGNPGLQAIWEDEKQRRRNRNESSQISQPESQDCRFVPATESEKQFQKRLQEVLKQNDFSVTLSGSVDYSNGSQEFSAELTLHSEILSPEMLPCSPANMIEVHKDTDLSKGNTKHKVEEALINEEAILNLIENSQTFQPLTQRLSETPVF/SHCAEMPTEVLPKKMTRPQKRKWSGPIIVCFLPISLYLS

MLEGKMADINFKEVTLIVSVVAACYWNSLFCGFVFDDVSAILDNKDLHPSTPLKTLFQNDFWGTPMSEERSHKSYRPLTVLTFRLNYLLSELKPMSYHLLNTVFHAVVSVIFLKVCRLFLDKRSSMIAALLFAVHPIHTEAVTGVVGRAELLSSVFFLAAFLSYTKSKGPDNSIVWTPIVLTVFLVAVATLCKEQGITVVGICCVYEVFVAQGYTLPMLCTVAGQFLRGKGSIPLSMLQTLVKLIVLMLSTLLLVVVRVQVIQSQLPVFTRFDNPAAVSPTPTRQLTFNYLLPVNAWLLLNPSELCCDWTMGTIPLIESFLDVRNLATFAFFCFLGALGIFSLRYPGDSSKTVLMALCLMALPFIPASNLFFPVGFVVAERVLYVPSMGFCILVAHGWQKISNKSVLKKLSWVCLSMVILTHALKTLHRNWDWESEYTLFMSALKVNKNNAKLWNNVGHALENEKNFEKALKYFLQATHVQPDDIGAHMNVGRTYKNLNRTREAEASYMLAKSLMPQIIPGKKYAARIAPNHLNVYINLANLIRANESRLEEADQLYRQAISMRPDFKQAYISRGELLLKMNKPLKAKEAYLKALELDRNNADLWYNLAIVYIELKEPNEALKNFNRALELNPKHKLALFNSAILMQESGEVKLRPEARKRLLNYVNEEPQDANGYFNLGMLAMDDKKDSEAESWMKKAIKLQPDFRSALFNLALLYSQTAKELKALPILEELLKYYPDHTKGLILKGDILMNQKKDIPGAKKCFEKILEMDPSNVQGKHNLCVVYFEEKELLKAERCLVETLALAPHEEYIQRHLSIVRDRISSSGIVEQPLAPADKTPGTEEREEIPSEDVKEISSESRPPQILKTNNNRNSKSNKQSTENADQDAPHKTTKDIKEIEKKRVAALKRLEEIERILNGE
MLEGKMADINFKEVTLIVSVVAACYWNSLFCGFVFDDVSAILDNKDLHPSTPLKTLFQNDFWGTPMSEERSHKSYRPLTVLTFRLNYLLSELKPMSYHLLNTVFHAVVSVIFLKVCRLFLDKRSSMIAALLFAVHPIHTEAVTGVVGRAELLSSVFFLAAFLSYTKSKGPDNSIVWTPIVLTVFLVAVATLCKEQGITVVGICCVYEVFVAQGYTLPMLCTVAGQFLRGKGSIPLSMLQTLVKLIVLMLSTLLLVVVRVQVIQSQLPVFTRFDNPAAVSPTPTRQLTFNYLLPVNAWLLLNPSELCCDWTMGTIPLIESFLDVRNLATFAFFCFLGALGIFSLRYPGDSSKTVLMALCLMALPFIPASNLFFPVGFVVAERVLYVPSMGFCILVAHGWQKISNKSVLKKLSWVCLSMVILTHALKTLHRNWDWESEYTLFMSALKVNKNNAKLWNNVGHALENEKNFEKALKYFLQATHVQPDDIGAHMNVGRTYKNLNRTREAEASYMLAKSLMPQIIPGKKYAARIAPNHLNVYINLANLIRANESRLEEADQLYRQAISMRPDFKQAYISRGELLLKMNKPLKAKEAYLKALELDRNNADLWYNLAIVYIELKEPNEALKNFNRALELNPKHKLALFNSAILMQESGK FPENVSI

