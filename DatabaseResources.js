// DatabaseResources.js - Defines biological database resources, example files, and fetch operations

/**
 * Database resources for fetch operations
 * Each database has a name, description, base URL, and a list of sub-databases
 * Each sub-database defines supported formats and URLs for searching and fetching
 */
export const DatabaseResources = {
  ncbi: {
    name: "NCBI",
    description: "National Center for Biotechnology Information",
    baseUrl: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    databases: [
      {
        id: "nucleotide",
        name: "Nucleotide",
        description: "DNA and RNA sequences",
        format: "genbank",
        searchUrl: "esearch.fcgi?db=nucleotide&term=",
        fetchUrl: "efetch.fcgi?db=nucleotide&id=$ID&rettype=gb&retmode=text",
        exampleIds: ["NC_045512.2", "MN908947", "U49845", "J01415.2"],
        altFormats: [
          { name: "FASTA", param: "rettype=fasta&retmode=text" },
          { name: "GenBank", param: "rettype=gb&retmode=text" },
          { name: "EMBL", param: "rettype=embl&retmode=text" }
        ]
      },
      {
        id: "protein",
        name: "Protein",
        description: "Protein sequences",
        format: "fasta",
        searchUrl: "esearch.fcgi?db=protein&term=",
        fetchUrl: "efetch.fcgi?db=protein&id=$ID&rettype=fasta&retmode=text",
        exampleIds: ["YP_009724390.1", "P05067", "NP_000509.1"],
        altFormats: [
          { name: "FASTA", param: "rettype=fasta&retmode=text" },
          { name: "GenBank", param: "rettype=gp&retmode=text" }
        ]
      },
      {
        id: "gene",
        name: "Gene",
        description: "Gene information",
        format: "text",
        searchUrl: "esearch.fcgi?db=gene&term=",
        fetchUrl: "efetch.fcgi?db=gene&id=$ID&retmode=text",
        exampleIds: ["672", "5468", "4609"]
      },
      {
        id: "pubmed",
        name: "PubMed",
        description: "Biomedical literature",
        format: "text",
        searchUrl: "esearch.fcgi?db=pubmed&term=",
        fetchUrl: "efetch.fcgi?db=pubmed&id=$ID&rettype=abstract&retmode=text",
        exampleIds: ["33027211", "34995640"]
      }
    ]
  },
  
  ensembl: {
    name: "Ensembl",
    description: "Genome browser for vertebrate genomes",
    baseUrl: "https://rest.ensembl.org/",
    databases: [
      {
        id: "gene",
        name: "Gene",
        description: "Gene information",
        format: "fasta",
        searchUrl: "lookup/symbol/homo_sapiens/$TERM?expand=1",
        fetchUrl: "sequence/id/$ID?type=genomic",
        exampleIds: ["ENSG00000139618", "ENSG00000141510", "ENSG00000012048"],
        altFormats: [
          { name: "Genomic FASTA", param: "type=genomic" },
          { name: "cDNA FASTA", param: "type=cdna" },
          { name: "CDS FASTA", param: "type=cds" }
        ]
      },
      {
        id: "transcript",
        name: "Transcript",
        description: "Transcript sequences",
        format: "fasta",
        searchUrl: "lookup/id/$TERM?expand=1",
        fetchUrl: "sequence/id/$ID?type=cdna",
        exampleIds: ["ENST00000380152", "ENST00000407559"],
        altFormats: [
          { name: "cDNA FASTA", param: "type=cdna" },
          { name: "CDS FASTA", param: "type=cds" },
          { name: "Protein FASTA", param: "type=protein" }
        ]
      },
      {
        id: "variant",
        name: "Variant",
        description: "Genetic variation information",
        format: "json",
        searchUrl: "variation/human/$TERM?",
        fetchUrl: "variation/human/$ID?",
        exampleIds: ["rs699", "rs6025", "rs28897743"]
      }
    ]
  },
  
  uniprot: {
    name: "UniProt",
    description: "Universal Protein Resource",
    baseUrl: "https://rest.uniprot.org/uniprotkb/",
    databases: [
      {
        id: "protein",
        name: "Protein",
        description: "Protein sequences and annotations",
        format: "fasta",
        searchUrl: "search?query=$TERM",
        fetchUrl: "$ID.fasta",
        exampleIds: ["P05067", "P53_HUMAN", "P42858"],
        altFormats: [
          { name: "FASTA", param: ".fasta" },
          { name: "Text", param: ".txt" },
          { name: "XML", param: ".xml" }
        ]
      },
      {
        id: "proteingb",
        name: "Protein (GenBank)",
        description: "Protein in GenBank format",
        format: "genbank",
        searchUrl: "search?query=$TERM",
        fetchUrl: "$ID.txt",
        exampleIds: ["P05067", "P53_HUMAN", "P42858"]
      },
      {
        id: "uniref",
        name: "UniRef",
        description: "Clustered sets of protein sequences",
        format: "fasta",
        searchUrl: "https://rest.uniprot.org/uniref/search?query=$TERM",
        fetchUrl: "https://rest.uniprot.org/uniref/$ID.fasta",
        exampleIds: ["UniRef90_P05067", "UniRef50_P42858"]
      }
    ]
  },
  
  pdb: {
    name: "PDB",
    description: "Protein Data Bank",
    baseUrl: "https://files.rcsb.org/download/",
    databases: [
      {
        id: "structure",
        name: "Structure",
        description: "3D macromolecular structures",
        format: "pdb",
        searchUrl: "https://data.rcsb.org/search/primitive/entry/$TERM",
        fetchUrl: "$ID.pdb",
        exampleIds: ["1HHO", "4HHB", "7KDT", "6VXX", "6M0J"],
        altFormats: [
          { name: "PDB", param: ".pdb" },
          { name: "mmCIF", param: ".cif" },
          { name: "FASTA", param: ".fasta" }
        ]
      },
      {
        id: "structurecif",
        name: "Structure (mmCIF)",
        description: "3D structures in mmCIF format",
        format: "mmcif",
        searchUrl: "https://data.rcsb.org/search/primitive/entry/$TERM",
        fetchUrl: "$ID.cif",
        exampleIds: ["1HHO", "4HHB", "7KDT", "6VXX", "6M0J"]
      },
      {
        id: "biologicalassembly",
        name: "Biological Assembly",
        description: "Biologically relevant assembly",
        format: "pdb",
        searchUrl: "https://data.rcsb.org/search/primitive/entry/$TERM",
        fetchUrl: "$ID-assembly1.pdb",
        exampleIds: ["1HHO", "4HHB", "7KDT"]
      }
    ]
  },
  
  ebi: {
    name: "EBI",
    description: "European Bioinformatics Institute",
    baseUrl: "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensembl&format=embl&id=",
    databases: [
      {
        id: "ensembl",
        name: "Ensembl (EMBL)",
        description: "EMBL format from Ensembl",
        format: "embl",
        searchUrl: "https://www.ensembl.org/Multi/Search/Results?q=$TERM",
        fetchUrl: "$ID",
        exampleIds: ["ENSG00000139618", "ENSG00000141510"],
        altFormats: [
          { name: "EMBL", param: "format=embl" },
          { name: "GenBank", param: "format=genbank" },
          { name: "FASTA", param: "format=fasta" }
        ]
      }
    ]
  },
  
  pfam: {
    name: "Pfam",
    description: "Protein family database",
    baseUrl: "https://pfam.xfam.org/",
    databases: [
      {
        id: "family",
        name: "Protein Family",
        description: "Protein domain families",
        format: "stockholm",
        searchUrl: "search/keyword?query=$TERM",
        fetchUrl: "family/$ID/alignment/full",
        exampleIds: ["PF00042", "PF00059", "PF00069"]
      },
      {
        id: "protein",
        name: "Protein",
        description: "Protein domain annotations",
        format: "json",
        searchUrl: "protein?entry=$TERM&output=xml",
        fetchUrl: "protein/$ID?output=xml",
        exampleIds: ["P05067", "P42858"]
      }
    ]
  },
  
  rfam: {
    name: "Rfam",
    description: "RNA family database",
    baseUrl: "https://rfam.xfam.org/",
    databases: [
      {
        id: "family",
        name: "RNA Family",
        description: "RNA families",
        format: "stockholm",
        searchUrl: "search?q=$TERM",
        fetchUrl: "family/$ID/alignment",
        exampleIds: ["RF00001", "RF00005", "RF00177"]
      }
    ]
  },
  
  treebase: {
    name: "TreeBASE",
    description: "Database of phylogenetic trees",
    baseUrl: "https://treebase.org/treebase-web/search/",
    databases: [
      {
        id: "study",
        name: "Study",
        description: "Phylogenetic studies",
        format: "newick",
        searchUrl: "study/find?query=$TERM",
        fetchUrl: "phylows/study/$ID?format=nexml",
        exampleIds: ["S1925", "S2484"]
      },
      {
        id: "tree",
        name: "Tree",
        description: "Individual phylogenetic trees",
        format: "newick",
        searchUrl: "tree/find?query=$TERM",
        fetchUrl: "phylows/tree/$ID?format=newick",
        exampleIds: ["TR12345", "TR67890"]
      }
    ]
  }
};

/**
 * Map database URL patterns to formats
 */
export const DatabaseUrlToFormat = {
  'eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide': 'genbank',
  'eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein': 'fasta',
  'rest.ensembl.org/sequence/id': 'fasta',
  'rest.uniprot.org/uniprotkb/': 'fasta',
  'files.rcsb.org/download/': 'pdb',
  'www.ebi.ac.uk/Tools/dbfetch/dbfetch': 'embl'
};

/**
 * Example files for testing and demonstration
 */
export const ExampleFiles = {
  fasta: {
    name: 'example.fasta',
    description: 'Example FASTA file with multiple protein sequences',
    content: `>sp|P62942|FKBP1A_HUMAN Peptidyl-prolyl cis-trans isomerase FKBP1A OS=Homo sapiens OX=9606 GN=FKBP1A PE=1 SV=2
MGVQVETISPGDGRTFPKRGQTCVVHYTGMLEDGKKFDSSRDRNKPFKFMLGKQEVIRG
WEEGVAQMSVGQRAKLTISPDYAYGATGHPGIIPPHATLVFDVELLKLE
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNP
KVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHH
FGKEFTPPVQAAYQKVVAGVANALAHKYH
>sp|P00533|EGFR_HUMAN Epidermal growth factor receptor OS=Homo sapiens OX=9606 GN=EGFR PE=1 SV=2
MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCE
VVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYA
LAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMS
MDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQC
AAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKS
CLLGNEGLPTLNEKDSGNITNSLKQFSGVGGEFVCAIMNTTPLISGLVTICPFSCRLAS
RQYQILTNTVESIQNRMLIESDGDGSCRPLNCSLTLQMVGFSRQDHTFDPYPLRDTWTY
LNATNQTTCSDGDTEAHQCLTSGSGGHNCNTKEVQCVNLNCTGCALHCLHSLGNRTYEC
GACGHTYTAPCAPCPHRCPEGTINSMPPSGCSPDHGRFQSKFDNGNSCTSHGCLNGGVC
MDHCLKCFGFNQDPAKHCADHPSNGGCRPTELPNQSQLCPPCHSGCILSARDPYCIWSG
SAQQCFTTKGVCNGMGCHSLDDSPQCHTCAHNFTGPNCEFPQGVCSSPDCKASCPLNYD
CEWNDSACDQCKKCCAGRCPGKCNGPKPQAEDKSIHQGTRRPPFTKTYGLFGDPHCEFS
ATRKNSLKIEGENIEEGTVTENSIPSEPGQTTLFLKAAELGRGYLSAIASGLTVLVLFG
AIFLMLVAGIYYRMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFE
KKVREGNSELKLIVNREERVENEPTGHLVQVAAYNAAKFAASCKDEANVDISAGGFLQA
PKGEKIPVACFLQNFSSIYRRTQDQEDGSLEDQDDEEDTLGGDGAKQRMAPAGRGGPRQ
RALPSAVVIDGITLPLEIKEHNEKGGRSSPTKKLSQNTEEKPSLGDTQPVETARTTSLR
EPRGEPHSHSSPLEPRDTRPLSRLLAPPTTPPSTCAPASSSTGALTPPDLLSHFNEEEE
LAFDEPVYANLNQAPLPPRADSFPQPAFSSGSLGPTPPGTESHLVERGAACMVTAVTCP
PLPPIMDLEDFPEPGACPPEPLSEATPRSPSTSRTPLLSSLSATSNNSTVACIDRNNAE
ERLQSQLLRRMLCAQTPAMRPNSNSGPPRAVDGGLKPSCPLPLMRQCSNQALNNLASEV
QIQNNRHSNSGPPEGPHCPPTAPEVPEPVMEAKKGSQPQKGPPEGPLKRLKLPKPSLEP
AQRPKTPKTEGFQPAEESEAEAPPGEKLVINEPQAAGGAGSVPPTLPNHTSTPEGRYVE
PVQNGHVAEEEERLQRALQQDLAAVPPRPNSLDFLQFQSDAGGMSPWNASVAPSPPPPN
PAEPHFAMHPTQAAEEYVEKPKMCPQDSAPVGPLGQVHVRKRVAIKPFQVSLCTPRMKI
NNSFLNESLASREFEVPVPSTPPSSKAGVTPFLDCQEQGRLNSMGQQNGQLLSLPTSWK
TLSPRHTSAMANTLKLLQLSHPHLVQAPHAPVPVPAPVGAGQVLCPQAPTGAGPNGNAV
SEYVNASVAPLSPPRDGCFGLAKAGPSASKSQTPLLSDLPRAPTPAPPEDVASTVGVAE
EPPTKKPRKMNGDPTSDRKANLSVDEFSVLKTPSGKQPRNNVALAVSRLVASSSSSGDS
SDRLPKSTTPPPLPNRDQAQFLGHCSKPAEPPKQFVPVPEHQNLTFDSNGHSMEKRHII
EALLVTELDKIKPIAYRELLDDLAHNCHSNDEVLAEAKILTQAANNGVTPLWAAVSPLK
`
  },
  fastq: {
    name: 'example.fastq',
    description: 'Example FASTQ file with sequencing reads',
    content: `@M01967:76:000000000-CFDPJ:1:1101:16889:1444 1:N:0:GCCAAT
TTTCCTGGGAGTAGGACACCATGTGACCCAAGCATGTTGGCCCAATGTTGCATCATTGATTCATCGACACTTCGTTGGTGTCTTTACGCTTGAGATATGTTGATACAACTGTGTTTTTGCCAGTTGTGGTTTGCACACTTTCAGCACTGTACTTGATCAAGA
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGG
@M01967:76:000000000-CFDPJ:1:1101:6168:1450 1:N:0:GCCAAT
ACTAGGTCGGGTGGTTCGGATGCTCCGGCGGGCAGGTACATCATATCTGGTGAATCTACCCAGGTCCGATTGGACGATCTCGAGGGGACCGTAAGGGCACGCAAGCTCGCTTTATCGCGAATCAGCGTAACTCTCTCCAGACGATAGACGGTGGCTCCGCC
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGFFGGGGGGGGGGGGGGGGGGGGGGGCFEGGGGGGGGGGGGGGGGGGEGGGGGF
@M01967:76:000000000-CFDPJ:1:1101:16458:1465 1:N:0:GCCAAT
ATCGCAACCTCAACGCGGATCTGCCAACCCGTCAGCACGCCCAGTACGGCTACAATGCTGCGGTCAAGTTCTTCGCCGAACTGGAAAATGACGCCTACATTGTCGAGTACCGCCGTAAGCAGACTCAGGTGGGGGAGTACGACGCGCTTTGCGTCGAGGT
+
CCCCCGGEGGDGGGGGGGFGGGGGGFGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFFGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGCFGGGGGGGGGGGGG
@M01967:76:000000000-CFDPJ:1:1101:17550:1468 1:N:0:GCCAAT
CAGACCGCCGATATTCTTGGCCTGAACCCAGAGCGTCTTCACCTGCAGCGCCTGGATCTTGCGGTGCGAGAAGTAGCGCGCGCTGTCGGGCTCGACGATGGTTTCGTTGCGGGCGGCGACACCCATGAGCATGATGCGCGGCGCGAAGCCTTCGCGGCTG
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@M01967:76:000000000-CFDPJ:1:1101:16410:1494 1:N:0:GCCAAT
TCGGACTCGAACCGCTCACCATCGCGTCGTTAGCGCGACGGACTCTAACCATTGAGTCTATGAGCAAATCGTCTCCCGGAGCAGCCTCTCCCGGAGCAGCCTCTCCCGGAGCGGCCGCCCCTACTACCGCCACCACCGCGCCTACCTACTGCGACGAGTG
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGCFGGGGGGGGGGGGGGGGFGGGGGGGFFGGGGFGGGGGGGGGGGFGGGGGGGGGGCFCFFGGGGGGFGGGGFGGGGGGGGG
`
  },
  genbank: {
    name: 'example.gb',
    description: 'Example GenBank file with COVID-19 spike protein',
    content: `LOCUS       YP_009724390             1273 aa            linear   VRL 15-JAN-2020
DEFINITION  surface glycoprotein [Severe acute respiratory syndrome coronavirus
            2].
ACCESSION   YP_009724390
VERSION     YP_009724390.1
DBLINK      BioProject: PRJNA485481
DBSOURCE    REFSEQ: accession NC_045512.2
KEYWORDS    RefSeq.
SOURCE      Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
  ORGANISM  Severe acute respiratory syndrome coronavirus 2
            Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes;
            Nidovirales; Cornidovirineae; Coronaviridae; Orthocoronavirinae;
            Betacoronavirus; Sarbecovirus.
REFERENCE   1  (residues 1 to 1273)
  AUTHORS   Wu,F., Zhao,S., Yu,B., Chen,Y.M., Wang,W., Song,Z.G., Hu,Y.,
            Tao,Z.W., Tian,J.H., Pei,Y.Y., Yuan,M.L., Zhang,Y.L., Dai,F.H.,
            Liu,Y., Wang,Q.M., Zheng,J.J., Xu,L., Holmes,E.C. and Zhang,Y.Z.
  TITLE     A new coronavirus associated with human respiratory disease in China
  JOURNAL   Nature 579 (7798), 265-269 (2020)
REFERENCE   2  (residues 1 to 1273)
  CONSRTM   NCBI Genome Project
  TITLE     Direct Submission
  JOURNAL   Submitted (28-JAN-2020) National Center for Biotechnology
            Information, NIH, Bethesda, MD 20894, USA
REFERENCE   3  (residues 1 to 1273)
  AUTHORS   Wu,F., Zhao,S., Yu,B., Chen,Y.M., Wang,W., Song,Z.G., Hu,Y.,
            Tao,Z.W., Tian,J.H., Pei,Y.Y., Yuan,M.L., Zhang,Y.L., Dai,F.H.,
            Liu,Y., Wang,Q.M., Zheng,J.J., Xu,L., Holmes,E.C. and Zhang,Y.Z.
  TITLE     Direct Submission
  JOURNAL   Submitted (05-JAN-2020) Shanghai Public Health Clinical Center &
            School of Public Health, Fudan University, Shanghai, China
COMMENT     PROVISIONAL REFSEQ: This record has not yet been subject to final
            NCBI review. The reference sequence is identical to QHD43416.
            Annotation was added using homology to SARS coronavirus.
            The spike protein of SARS-CoV-2 contains a receptor-binding domain
            (RBD) that interacts with the angiotensin-converting enzyme 2 (ACE2)
            on host cells. For more information on the spike proteins of
            sarbecoviruses and the differences between SARS-CoV-1 and
            SARS-CoV-2, see PMIDs 35817756 and 32668444.
            
            ##Assembly-Data-START##
            Assembly Method       :: megahit v. V1.1.3
            Sequencing Technology :: Illumina
            ##Assembly-Data-END##
FEATURES             Location/Qualifiers
     source          1..1273
                     /organism="Severe acute respiratory syndrome coronavirus 2"
                     /isolate="Wuhan-Hu-1"
                     /host="Homo sapiens"
                     /db_xref="taxon:2697049"
                     /country="China"
                     /collection_date="Dec-2019"
     Protein         1..1273
                     /product="surface glycoprotein"
                     /note="S glycoprotein; structural protein; spike protein"
     Region          1..13
                     /region_name="Signal"
                     /note="propagated from UniProtKB/Swiss-Prot (P0DTC2.1)"
                     /db_xref="CDD:340937"
     Region          14..1273
                     /region_name="Mature chain"
                     /note="propagated from UniProtKB/Swiss-Prot (P0DTC2.1)"
                     /db_xref="CDD:340937"
     Region          319..541
                     /region_name="receptor binding domain"
                     /note="propagated from UniProtKB/Swiss-Prot (P0DTC2.1)"
                     /db_xref="CDD:340937"
     Site            order(493,501..502,505)
                     /site_type="other"
                     /note="ACE2 binding site [polypeptide binding]"
                     /db_xref="CDD:340937"
     Site            order(417,446..447,449,453,455..456,475..476,478..489,
                     493..494,496..498,500..505,516)
                     /site_type="other"
                     /note="receptor binding motif [polypeptide binding]"
                     /db_xref="CDD:340937"
     Site            order(417,421..424,455,459..461,463,469..472,479,486..491,
                     493..494,496,498,500..501,503..505)
                     /site_type="other"
                     /note="RBD interface [polypeptide binding]"
                     /db_xref="CDD:340937"
     Site            682..685
                     /site_type="other"
                     /note="furin cleavage site"
                     /db_xref="CDD:340937"
     Site            order(816,896,943,970)
                     /site_type="other"
                     /note="N-glycosylation site"
                     /db_xref="CDD:340937"
     Site            816
                     /site_type="other"
                     /note="N-linked glycosylation site"
                     /db_xref="CDD:340937"
     Region          818..1197
                     /region_name="Corona_S2"
                     /note="Corona S2 glycoprotein; pfam01601"
                     /db_xref="CDD:279881"
     Region          910..988
                     /region_name="heptad repeat 1"
                     /note="propagated from UniProtKB/Swiss-Prot (P0DTC2.1)"
                     /db_xref="CDD:340937"
     Region          1162..1203
                     /region_name="heptad repeat 2"
                     /note="propagated from UniProtKB/Swiss-Prot (P0DTC2.1)"
                     /db_xref="CDD:340937"
     Site            1214
                     /site_type="other"
                     /note="Palmitoylated Cys [lipid moiety-binding]"
                     /db_xref="CDD:340937"
     Site            order(1226,1235..1236,1244,1247,1248..1250,1253,1261)
                     /site_type="other"
                     /note="Palmitoylated Cys [lipid moiety-binding]"
                     /db_xref="CDD:340937"
     Site            order(1114,1118,1134,1158,1173..1175,1182,1191,1195,1199,
                     1235,1236,1240,1243,1247,1254,1264)
                     /site_type="other"
                     /note="S2' interface [polypeptide binding]"
                     /db_xref="CDD:340937"
ORIGIN      
        1 mfvflvllpl vssqcvnltt rtqlppaytn sftrgvyypd kvfrssvlhs tqdlflpffs
       61 nvtwfhaihv sgtngtkrfd npvlpfndgv yfasteksni irgwifgttl dsktqslliv
      121 nnatnvvikv cefqfcndpf lgvyyhknnk swmesefrvy ssannctfey vsqpflmdle
      181 gkqgnfknlr efvfknidgy fkiyskhtpi nlvrdlpqgf saleplvdlp iginitrfqt
      241 llalhrsylt pgdsssgwta gaaayyvgyl qprtfllkyn engtitdavd caldplsetk
      301 ctlksftvek giyqtsnfrv qptesivrfp nitnlcpfge vfnatrfasv yawnrkrisn
      361 cvadysvlyn sasfstfkcy gvsptklndl cftnvyadsf virgdevrqi apgqtgkiad
      421 ynyklpddft gcviawnsnn ldskvggnyn ylyrlfrksnlkpferdist eiyqagstpc
      481 ngvegfncyf plqsygfqpt ngvgyqpyrv vvlsfellha patvcgpkks tnlvknkcvn
      541 fnfngltgtg vltesnkkfl pfqqfgrdia dttdavrdpq tleilditpc sfggvsvitp
      601 gtntsnqvav lyqdvnctev pvaihadqlt ptwrvystgs nvfqtragcl igaehvnnsy
      661 ecdipigagi casyqtqtns prrarsvasq siiaytmslg aensvaysnn siaiptnfti
      721 svtteilpvs mtktsvdctm yicgdstecs nlllqygsfc tqlnraltgi aveqdkntqe
      781 vfaqvkqiyk tppikdfggf nfsqilpdps kpskrsfied llfnkvtlad agfikqygdc
      841 lgdiaardli caqkfngltv lpplltdemi aqytsallag titsgwtfga gaalqipfam
      901 qmayrfngig vtqnvlyenq klianqfnsa igkiqdslss tasalgklqd vvnqnaqaln
      961 tlvkqlssnf gaissvlndi lsrldkveae vqidrlitgr lqslqtyvtq qliraaeira
     1021 sanlaatkms ecvlgqskrv dfcgkgyhlm sfpqsaphgv vflhvtyvpa qeknfttapa
     1081 ichdgkahfp regvfvsngt hwfvtqrnfy epqiittdnt fvsgncdvvi givnntvydp
     1141 lqpeldsfke eldkyfknht spdvdlgdis ginasvvniq keidrlneva knlneslidl
     1201 qelgkyeqyi kwpwyiwlgf iagliaivmv timlccmtsc csclkgccsc gscckfdedd
     1261 sepvlkgvkl hyt
//`
  },
  pdb: {
    name: 'example.pdb',
    description: 'Example PDB file with hemoglobin alpha chain',
    content: `HEADER    OXYGEN TRANSPORT                        14-JUL-83   1HHO      
TITLE     THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 ANGSTROMS     
TITLE    2 RESOLUTION                                                           
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: HEMOGLOBIN ALPHA CHAIN;                                    
COMPND   3 CHAIN: A, C;                                                         
COMPND   4 ENGINEERED: YES                                                      
COMPND   5 MOL_ID: 2;                                                           
COMPND   6 MOLECULE: HEMOGLOBIN BETA CHAIN;                                     
COMPND   7 CHAIN: B, D;                                                         
COMPND   8 ENGINEERED: YES                                                      
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;                                   
SOURCE   3 ORGANISM_COMMON: HUMAN;                                              
SOURCE   4 ORGANISM_TAXID: 9606;                                                
SOURCE   5 GENE: HBA;                                                           
SOURCE   6 TISSUE: BLOOD                                                        
SOURCE   7 MOL_ID: 2;                                                           
SOURCE   8 ORGANISM_SCIENTIFIC: HOMO SAPIENS;                                   
SOURCE   9 ORGANISM_COMMON: HUMAN;                                              
SOURCE  10 ORGANISM_TAXID: 9606;                                                
SOURCE  11 GENE: HBB;                                                           
SOURCE  12 TISSUE: BLOOD                                                        
KEYWDS    OXYGEN TRANSPORT, HAEM, OXYGEN BINDING PROTEIN                         
EXPDTA    X-RAY DIFFRACTION                                                     
AUTHOR    G.FERMI,M.F.PERUTZ                                                    
HELIX    1   1 VAL A   1  LEU A  18  1                                  18    
HELIX    2   2 PRO A  36  THR A  41  1                                   6    
HELIX    3   3 PHE A  43  ASP A  54  1                                  12    
HELIX    4   4 LYS A  60  LYS A  61  1                                   2    
HELIX    5   5 LEU A  76  GLU A  85  1                                  10    
HELIX    6   6 LYS A  90  LEU A  97  1                                   8    
HELIX    7   7 THR A 102  HIS A 122  1                                  21    
HELIX    8   8 ALA A 130  ALA A 137  1                                   8    
HELIX    9   1 HIS B   2  SER B  17  1                                  16    
ATOM      1  N   VAL A   1      11.804  18.255  16.077  1.00 19.49           N  
ATOM      2  CA  VAL A   1      11.629  19.074  14.856  1.00 18.25           C  
ATOM      3  C   VAL A   1      10.477  18.607  13.991  1.00 14.92           C  
ATOM      4  O   VAL A   1       9.467  18.112  14.482  1.00 16.47           O  
ATOM      5  CB  VAL A   1      11.628  20.577  15.152  1.00 22.78           C  
ATOM      6  CG1 VAL A   1      11.501  21.362  13.861  1.00 24.66           C  
ATOM      7  CG2 VAL A   1      12.895  20.966  15.895  1.00 24.39           C  
ATOM      8  N   LEU A   2      10.647  18.742  12.686  1.00 14.16           N  
ATOM      9  CA  LEU A   2       9.679  18.314  11.678  1.00 13.20           C  
ATOM     10  C   LEU A   2       9.370  19.468  10.715  1.00 12.25           C  
ATOM     11  O   LEU A   2      10.111  20.442  10.616  1.00 11.93           O  
ATOM     12  CB  LEU A   2      10.184  17.096  10.889  1.00 15.94           C  
ATOM     13  CG  LEU A   2      11.328  17.326   9.888  1.00 22.66           C  
ATOM     14  CD1 LEU A   2      12.659  17.354  10.620  1.00 25.00           C  
ATOM     15  CD2 LEU A   2      11.248  16.283   8.801  1.00 23.84           C  
ATOM     16  N   SER A   3       8.262  19.355  10.008  1.00 10.04           N  
ATOM     17  CA  SER A   3       7.868  20.383   9.046  1.00 10.38           C  
ATOM     18  C   SER A   3       8.764  20.361   7.824  1.00 11.92           C  
ATOM     19  O   SER A   3       9.145  19.286   7.371  1.00 14.16           O  
ATOM     20  CB  SER A   3       6.416  20.162   8.637  1.00 10.71           C  
ATOM     21  OG  SER A   3       6.339  19.131   7.671  1.00 14.08           O  
ATOM     22  N   PRO A   4       9.106  21.533   7.283  1.00 10.45           N  
ATOM     23  CA  PRO A   4      10.039  21.667   6.163  1.00 10.20           C  
ATOM     24  C   PRO A   4       9.672  20.925   4.883  1.00 11.23           C  
ATOM     25  O   PRO A   4      10.511  20.303   4.241  1.00 12.70           O  
ATOM     26  CB  PRO A   4      10.103  23.182   5.937  1.00 12.14           C  
TER    2177      ILE D 146                                                      
HETATM 2178  FE  HEM A 142      10.846  16.960  20.168  1.00 12.18          FE  
HETATM 2179  CHA HEM A 142      12.242  17.313  17.340  1.00 10.37           C  
HETATM 2180  CHB HEM A 142       8.329  15.958  19.344  1.00 10.63           C  
HETATM 2181  CHC HEM A 142       9.476  16.728  23.004  1.00 11.40           C  
HETATM 2182  CHD HEM A 142      13.379  18.018  21.043  1.00 13.32           C  
HETATM 2183  NA  HEM A 142      10.587  16.508  18.239  1.00 12.49           N  
HETATM 2184  C1A HEM A 142      11.501  16.918  17.275  1.00 11.70           C  
HETATM 2185  C2A HEM A 142      11.213  16.823  15.878  1.00 10.93           C  
HETATM 2186  C3A HEM A 142      10.075  16.259  15.833  1.00 14.80           C  
HETATM 2187  C4A HEM A 142       9.593  16.091  17.207  1.00 12.16           C  
HETATM 2188  CMA HEM A 142       9.300  15.814  14.622  1.00 13.08           C  
HETATM 2189  CAA HEM A 142      12.043  17.118  14.667  1.00 15.37           C  
HETATM 2190  CBA HEM A 142      13.316  17.874  14.987  1.00 17.36           C  
HETATM 2191  CGA HEM A 142      14.144  18.300  13.802  1.00 22.14           C  
HETATM 2192  O1A HEM A 142      13.951  17.804  12.682  1.00 25.33           O  
HETATM 2193  O2A HEM A 142      15.078  19.056  13.983  1.00 24.74           O  
HETATM 2194  NB  HEM A 142       9.091  16.317  20.022  1.00 10.45           N  
CONECT 2178 2183 2194 2212 2224                                                 
CONECT 2179 2184 2266                                                           
CONECT 2180 2187 2195                                                           
CONECT 2181 2202 2208                                                           
CONECT 2182 2211 2225                                                           
CONECT 2183 2178 2184 2187                                                      
CONECT 2184 2179 2183 2185                                                      
CONECT 2185 2184 2186 2189                                                      
CONECT 2186 2185 2187 2188                                                      
CONECT 2187 2180 2183 2186                                                      
CONECT 2188 2186                                                                
CONECT 2189 2185 2190                                                           
CONECT 2190 2189 2191                                                           
CONECT 2191 2190 2192 2193                                                      
CONECT 2192 2191                                                                
CONECT 2193 2191                                                                
CONECT 2194 2178 2195 2202                                                      
MASTER      332    0    1    1    5    0    0    6 2238    4   40   21          
END                                                                             
`
  },
  clustal: {
    name: 'example.clustal',
    description: 'Example CLUSTAL alignment of globin proteins',
    content: `CLUSTAL W (1.83) multiple sequence alignment

MYG_PHYCA      ---------VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKT
HBA_HUMAN      VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK
HBB_HUMAN      VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
GLB5_PETMA     PIVDTGSVAPLSAAEKTKIRSAWAPVYSTYETSGVDILVKFFTSTPAAQEFFPKFKGLTT
                         .    .  . .     ..    . .        .  .  . .  .    . 

MYG_PHYCA      EAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAI
HBA_HUMAN      KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA
HBB_HUMAN      KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK
GLB5_PETMA     ADQLKKSADVRWHAERIINAVNDAVASMDDTEKMSMKLRDLSGKHAKSFQVDPQYFKVLA
                    .. .   . .  . . .   .. .   .    .   .  . . .    . .   . 

MYG_PHYCA      IHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG
HBA_HUMAN      VHASLDKFLASVSTVLTSKYR----------
HBB_HUMAN      EFTPPVQAAYQKVVAGVANALAHKYH------
GLB5_PETMA     AVI-ADTVAAGDAGFEKLMSMICILLRSAY--
                . ..       .  . .   .          `
  }
};

/**
 * Generate a URL for downloading a sequence from a database
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @param {string} accessionId - The accession ID
 * @returns {string} Generated URL
 */
export const generateDatabaseUrl = (databaseId, subDatabaseId, accessionId) => {
  const database = DatabaseResources[databaseId];
  if (!database) return null;
  
  const subDatabase = database.databases.find(db => db.id === subDatabaseId);
  if (!subDatabase) return null;
  
  return database.baseUrl + subDatabase.fetchUrl.replace('$ID', accessionId);
};

/**
 * Search for a term in a database
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @param {string} searchTerm - The search term
 * @returns {string} Search URL
 */
export const generateSearchUrl = (databaseId, subDatabaseId, searchTerm) => {
  const database = DatabaseResources[databaseId];
  if (!database) return null;
  
  const subDatabase = database.databases.find(db => db.id === subDatabaseId);
  if (!subDatabase) return null;
  
  // For NCBI, construct a complete URL including the base URL
  if (databaseId === 'ncbi') {
    return database.baseUrl + subDatabase.searchUrl + encodeURIComponent(searchTerm);
  } 
  // For Ensembl, replace $TERM with the search term
  else if (databaseId === 'ensembl') {
    return database.baseUrl + subDatabase.searchUrl.replace('$TERM', encodeURIComponent(searchTerm));
  }
  // For other databases
  else {
    return database.baseUrl + subDatabase.searchUrl.replace('$TERM', encodeURIComponent(searchTerm));
  }
};

/**
 * Detect if a URL is from a supported database and extract the format
 * @param {string} url - The URL to analyze
 * @returns {string} The detected format or null
 */
export const detectDatabaseFormat = (url) => {
  if (!url) return null;
  
  // Check each pattern
  for (const [pattern, format] of Object.entries(DatabaseUrlToFormat)) {
    if (url.includes(pattern)) {
      return format;
    }
  }
  
  // Check file extensions
  if (url.endsWith('.fasta') || url.endsWith('.fa')) {
    return 'fasta';
  } else if (url.endsWith('.gb') || url.endsWith('.gbk')) {
    return 'genbank';
  } else if (url.endsWith('.embl')) {
    return 'embl';
  } else if (url.endsWith('.pdb')) {
    return 'pdb';
  } else if (url.endsWith('.cif')) {
    return 'mmcif';
  }
  
  return null;
};

/**
 * Extract accession ID from a database URL
 * @param {string} url - The URL to analyze
 * @returns {string} The extracted accession ID or null
 */
export const extractAccessionId = (url) => {
  if (!url) return null;
  
  // Extract from NCBI URL
  const ncbiMatch = url.match(/id=([^&]+)/);
  if (ncbiMatch) return ncbiMatch[1];
  
  // Extract from Ensembl URL
  const ensemblMatch = url.match(/id\/([^?\/]+)/);
  if (ensemblMatch) return ensemblMatch[1];
  
  // Extract from UniProt URL
  const uniprotMatch = url.match(/uniprot\.org\/uniprotkb\/([^.]+)/);
  if (uniprotMatch) return uniprotMatch[1];
  
  // Extract from PDB URL
  const pdbMatch = url.match(/download\/([^.]+)/);
  if (pdbMatch) return pdbMatch[1];
  
  // Extract from URL path
  const pathMatch = url.match(/\/([A-Za-z0-9_]+)\.(fasta|gb|pdb|embl|cif)$/);
  if (pathMatch) return pathMatch[1];
  
  return null;
};

/**
 * Fetch data from a database URL
 * @param {string} url - The URL to fetch from
 * @param {function} successCallback - Callback for successful fetch
 * @param {function} errorCallback - Callback for error handling
 */
export const fetchFromUrl = async (url, successCallback, errorCallback) => {
  try {
    // Extract database information
    const { database, subDatabase } = getDatabaseFromUrl(url);
    
    // Determine expected format
    let format = null;
    if (database && subDatabase) {
      const dbResource = DatabaseResources[database];
      const subDbResource = dbResource.databases.find(db => db.id === subDatabase);
      if (subDbResource) {
        format = subDbResource.format;
      }
    }
    
    // If not determined from database info, try to detect from URL
    if (!format) {
      format = detectDatabaseFormat(url);
    }
    
    // Fetch the data
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`Failed to fetch from ${url}: ${response.status} ${response.statusText}`);
    }
    
    const content = await response.text();
    
    // Call success callback with the content and detected format
    successCallback({
      content,
      format,
      url,
      database,
      subDatabase,
      accessionId: extractAccessionId(url)
    });
    
  } catch (error) {
    errorCallback(`Error fetching data: ${error.message}`);
  }
};

/**
 * Search a database for a term
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @param {string} searchTerm - The search term
 * @param {function} successCallback - Callback for successful search
 * @param {function} errorCallback - Callback for error handling
 */
export const searchDatabase = async (databaseId, subDatabaseId, searchTerm, successCallback, errorCallback) => {
  try {
    // Generate search URL
    const searchUrl = generateSearchUrl(databaseId, subDatabaseId, searchTerm);
    if (!searchUrl) {
      throw new Error('Could not generate search URL');
    }
    
    // Fetch search results
    const response = await fetch(searchUrl);
    if (!response.ok) {
      throw new Error(`Search failed: ${response.status} ${response.statusText}`);
    }
    
    // Parse response based on database
    let results = [];
    const contentType = response.headers.get('Content-Type') || '';
    const responseText = await response.text();
    
    if (databaseId === 'ncbi') {
      results = parseNcbiSearchResults(responseText, subDatabaseId);
    } else if (databaseId === 'ensembl') {
      results = parseEnsemblSearchResults(responseText, subDatabaseId);
    } else if (databaseId === 'uniprot') {
      results = parseUniprotSearchResults(responseText, subDatabaseId);
    } else if (databaseId === 'pdb') {
      results = parsePdbSearchResults(responseText, subDatabaseId);
    } else {
      // Generic parsing based on content type
      if (contentType.includes('json')) {
        try {
          const jsonData = JSON.parse(responseText);
          results = parseGenericJsonResults(jsonData, databaseId, subDatabaseId);
        } catch (e) {
          console.error('Error parsing JSON:', e);
        }
      } else {
        // Try to extract IDs from HTML or text response
        results = extractIdsFromText(responseText, databaseId, subDatabaseId);
      }
    }
    
    // Call success callback with the results
    successCallback({
      term: searchTerm,
      database: databaseId,
      subDatabase: subDatabaseId,
      results,
      raw: responseText
    });
    
  } catch (error) {
    errorCallback(`Search error: ${error.message}`);
  }
};

/**
 * Parse NCBI search results
 * @param {string} responseText - The search response
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of search result objects
 */
const parseNcbiSearchResults = (responseText, subDatabaseId) => {
  const results = [];
  
  // NCBI eSearch returns XML
  const idMatch = responseText.match(/<Id>(\d+)<\/Id>/g);
  if (idMatch) {
    const ids = idMatch.map(id => id.replace(/<Id>|<\/Id>/g, ''));
    
    // Create a result object for each ID
    ids.forEach(id => {
      results.push({
        id,
        title: `${subDatabaseId.toUpperCase()} entry ${id}`,
        database: 'ncbi',
        subDatabase: subDatabaseId
      });
    });
  }
  
  return results;
};

/**
 * Parse Ensembl search results
 * @param {string} responseText - The search response
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of search result objects
 */
const parseEnsemblSearchResults = (responseText, subDatabaseId) => {
  const results = [];
  
  try {
    // Ensembl REST API returns JSON
    const jsonData = JSON.parse(responseText);
    
    if (jsonData.id) {
      // Single result
      results.push({
        id: jsonData.id,
        title: jsonData.display_name || jsonData.id,
        description: jsonData.description || '',
        database: 'ensembl',
        subDatabase: subDatabaseId
      });
    } else if (Array.isArray(jsonData)) {
      // Array of results
      jsonData.forEach(item => {
        results.push({
          id: item.id,
          title: item.display_name || item.id,
          description: item.description || '',
          database: 'ensembl',
          subDatabase: subDatabaseId
        });
      });
    }
  } catch (e) {
    console.error('Error parsing Ensembl results:', e);
  }
  
  return results;
};

/**
 * Parse UniProt search results
 * @param {string} responseText - The search response
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of search result objects
 */
const parseUniprotSearchResults = (responseText, subDatabaseId) => {
  const results = [];
  
  try {
    // UniProt REST API returns JSON
    const jsonData = JSON.parse(responseText);
    
    if (jsonData.results && Array.isArray(jsonData.results)) {
      jsonData.results.forEach(item => {
        results.push({
          id: item.primaryAccession,
          title: item.uniProtkbId || item.primaryAccession,
          description: item.proteinDescription?.recommendedName?.fullName?.value || '',
          organism: item.organism?.scientificName || '',
          database: 'uniprot',
          subDatabase: subDatabaseId
        });
      });
    }
  } catch (e) {
    console.error('Error parsing UniProt results:', e);
  }
  
  return results;
};

/**
 * Parse PDB search results
 * @param {string} responseText - The search response
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of search result objects
 */
const parsePdbSearchResults = (responseText, subDatabaseId) => {
  const results = [];
  
  try {
    // PDB API returns JSON
    const jsonData = JSON.parse(responseText);
    
    if (jsonData.result_set && Array.isArray(jsonData.result_set)) {
      jsonData.result_set.forEach(item => {
        results.push({
          id: item.identifier,
          title: item.title || item.identifier,
          description: item.description || '',
          database: 'pdb',
          subDatabase: subDatabaseId
        });
      });
    }
  } catch (e) {
    console.error('Error parsing PDB results:', e);
    
    // Fallback: Try to extract PDB IDs from text
    const pdbIdRegex = /\b[1-9][A-Za-z0-9]{3}\b/g;
    const matches = responseText.match(pdbIdRegex);
    
    if (matches) {
      const uniqueIds = [...new Set(matches)];
      uniqueIds.forEach(id => {
        results.push({
          id,
          title: `PDB structure ${id}`,
          database: 'pdb',
          subDatabase: subDatabaseId
        });
      });
    }
  }
  
  return results;
};

/**
 * Parse generic JSON search results
 * @param {object} jsonData - The parsed JSON data
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of search result objects
 */
const parseGenericJsonResults = (jsonData, databaseId, subDatabaseId) => {
  const results = [];
  
  // Try to extract IDs and titles from common JSON structures
  if (Array.isArray(jsonData)) {
    jsonData.forEach(item => {
      const id = item.id || item.accession || item.identifier || '';
      if (id) {
        results.push({
          id,
          title: item.name || item.title || item.description || id,
          description: item.description || '',
          database: databaseId,
          subDatabase: subDatabaseId
        });
      }
    });
  } else if (jsonData.results && Array.isArray(jsonData.results)) {
    jsonData.results.forEach(item => {
      const id = item.id || item.accession || item.identifier || '';
      if (id) {
        results.push({
          id,
          title: item.name || item.title || item.description || id,
          description: item.description || '',
          database: databaseId,
          subDatabase: subDatabaseId
        });
      }
    });
  } else if (jsonData.hits && Array.isArray(jsonData.hits)) {
    jsonData.hits.forEach(item => {
      const id = item.id || item.accession || item.identifier || '';
      if (id) {
        results.push({
          id,
          title: item.name || item.title || item.description || id,
          description: item.description || '',
          database: databaseId,
          subDatabase: subDatabaseId
        });
      }
    });
  }
  
  return results;
};

/**
 * Extract IDs from text response
 * @param {string} text - The text to extract IDs from
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of search result objects
 */
const extractIdsFromText = (text, databaseId, subDatabaseId) => {
  const results = [];
  
  // Define regex patterns for different databases
  const patterns = {
    'ncbi': {
      'nucleotide': /\b[A-Z]{1,2}\_\d+(\.\d+)?\b/g,
      'protein': /\b[A-Z]{2}\_\d+(\.\d+)?\b/g,
      'gene': /\b\d{1,7}\b/g,
      'pubmed': /\b\d{5,9}\b/g
    },
    'ensembl': {
      'gene': /\bENS[A-Z]*G\d{11}\b/g,
      'transcript': /\bENS[A-Z]*T\d{11}\b/g,
      'variant': /\brs\d+\b/g
    },
    'uniprot': {
      'protein': /\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b/g,
      'proteingb': /\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b/g,
      'uniref': /\bUniRef\d+\_[OPQ][0-9][A-Z0-9]{3}[0-9]|UniRef\d+\_[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b/g
    },
    'pdb': {
      'structure': /\b[1-9][A-Za-z0-9]{3}\b/g,
      'structurecif': /\b[1-9][A-Za-z0-9]{3}\b/g,
      'biologicalassembly': /\b[1-9][A-Za-z0-9]{3}\b/g
    },
    'pfam': {
      'family': /\bPF\d{5}\b/g,
      'protein': /\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b/g
    },
    'rfam': {
      'family': /\bRF\d{5}\b/g
    },
    'treebase': {
      'study': /\bS\d+\b/g,
      'tree': /\bTR\d+\b/g
    }
  };
  
  // Get the appropriate regex pattern
  const pattern = patterns[databaseId]?.[subDatabaseId] || /\b[A-Za-z0-9_]+\b/g;
  
  // Extract matches
  const matches = text.match(pattern);
  
  if (matches) {
    // Remove duplicates
    const uniqueIds = [...new Set(matches)];
    
    // Create result objects
    uniqueIds.forEach(id => {
      results.push({
        id,
        title: `${databaseId.toUpperCase()} ${subDatabaseId} entry: ${id}`,
        database: databaseId,
        subDatabase: subDatabaseId
      });
    });
  }
  
  return results;
};

/**
 * Generate a URL with alternate format
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @param {string} accessionId - The accession ID
 * @param {string} altFormat - The alternate format parameter
 * @returns {string} Generated URL
 */
export const generateAltFormatUrl = (databaseId, subDatabaseId, accessionId, altFormat) => {
  const database = DatabaseResources[databaseId];
  if (!database) return null;
  
  const subDatabase = database.databases.find(db => db.id === subDatabaseId);
  if (!subDatabase) return null;
  
  const altFormatObj = subDatabase.altFormats?.find(fmt => fmt.param === altFormat);
  if (!altFormatObj) return null;
  
  // Start with base URL
  let url = database.baseUrl;
  
  // For NCBI, replace the parameters
  if (databaseId === 'ncbi') {
    const baseUrl = url + subDatabase.fetchUrl.split('&rettype=')[0] + '&';
    return baseUrl + altFormat;
  } 
  // For Ensembl, replace the type parameter
  else if (databaseId === 'ensembl') {
    const fetchUrl = subDatabase.fetchUrl.replace('$ID', accessionId);
    return database.baseUrl + fetchUrl.split('?')[0] + '?' + altFormat;
  }
  // For UniProt, replace the extension
  else if (databaseId === 'uniprot') {
    if (altFormat.startsWith('.')) {
      return database.baseUrl + accessionId + altFormat;
    } else {
      return database.baseUrl + accessionId + '?' + altFormat;
    }
  }
  // For PDB, replace the extension
  else if (databaseId === 'pdb') {
    if (altFormat.startsWith('.')) {
      return database.baseUrl + accessionId + altFormat;
    } else {
      return database.baseUrl + accessionId + '?' + altFormat;
    }
  }
  // For EBI, replace the format parameter
  else if (databaseId === 'ebi') {
    const baseUrl = database.baseUrl.replace('format=embl', '');
    return baseUrl + altFormat + '&id=' + accessionId;
  }
  
  return null;
};

/**
 * Get alternative format options for a database
 * @param {string} databaseId - The database ID
 * @param {string} subDatabaseId - The sub-database ID
 * @returns {Array} Array of alternative format objects
 */
export const getAltFormatOptions = (databaseId, subDatabaseId) => {
  const database = DatabaseResources[databaseId];
  if (!database) return [];
  
  const subDatabase = database.databases.find(db => db.id === subDatabaseId);
  if (!subDatabase || !subDatabase.altFormats) return [];
  
  return subDatabase.altFormats;
};

/**
 * Get database and sub-database from URL
 * @param {string} url - The URL to analyze
 * @returns {object} Object with database and subDatabase IDs
 */
export const getDatabaseFromUrl = (url) => {
  if (!url) return { database: null, subDatabase: null };
  
  if (url.includes('ncbi.nlm.nih.gov')) {
    const dbMatch = url.match(/db=([^&]+)/);
    if (dbMatch) {
      const dbName = dbMatch[1];
      if (dbName === 'nucleotide') return { database: 'ncbi', subDatabase: 'nucleotide' };
      if (dbName === 'protein') return { database: 'ncbi', subDatabase: 'protein' };
      if (dbName === 'gene') return { database: 'ncbi', subDatabase: 'gene' };
      if (dbName === 'pubmed') return { database: 'ncbi', subDatabase: 'pubmed' };
    }
  } else if (url.includes('rest.ensembl.org')) {
    if (url.includes('/sequence/id/')) {
      if (url.includes('type=cdna')) return { database: 'ensembl', subDatabase: 'transcript' };
      return { database: 'ensembl', subDatabase: 'gene' };
    }
    if (url.includes('/variation/human/')) {
      return { database: 'ensembl', subDatabase: 'variant' };
    }
  } else if (url.includes('uniprot.org')) {
    if (url.includes('uniref')) {
      return { database: 'uniprot', subDatabase: 'uniref' };
    } else {
      if (url.endsWith('.fasta')) return { database: 'uniprot', subDatabase: 'protein' };
      if (url.endsWith('.txt')) return { database: 'uniprot', subDatabase: 'proteingb' };
      return { database: 'uniprot', subDatabase: 'protein' };
    }
  } else if (url.includes('rcsb.org') || url.includes('files.rcsb.org')) {
    if (url.includes('-assembly')) return { database: 'pdb', subDatabase: 'biologicalassembly' };
    else if (url.endsWith('.pdb')) return { database: 'pdb', subDatabase: 'structure' };
    else if (url.endsWith('.cif')) return { database: 'pdb', subDatabase: 'structurecif' };
  } else if (url.includes('ebi.ac.uk')) {
    if (url.includes('dbfetch')) return { database: 'ebi', subDatabase: 'ensembl' };
  } else if (url.includes('pfam.xfam.org')) {
    if (url.includes('/family/')) return { database: 'pfam', subDatabase: 'family' };
    if (url.includes('/protein/')) return { database: 'pfam', subDatabase: 'protein' };
  } else if (url.includes('rfam.xfam.org')) {
    return { database: 'rfam', subDatabase: 'family' };
  } else if (url.includes('treebase.org')) {
    if (url.includes('/study/')) return { database: 'treebase', subDatabase: 'study' };
    if (url.includes('/tree/')) return { database: 'treebase', subDatabase: 'tree' };
  }
  
  return { database: null, subDatabase: null };
};
