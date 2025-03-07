// FormatConverters.js - Handles conversions between different biological data formats

/**
 * Main function to convert a file from one format to another
 * @param {string} content - The content of the file to convert
 * @param {string} sourceFormat - The source format key
 * @param {string} targetFormat - The target format key
 * @param {object} options - Conversion options
 * @param {function} progressCallback - Callback for progress updates
 * @param {function} successCallback - Callback for successful conversion
 * @param {function} errorCallback - Callback for error handling
 */
export const convertFile = (
  content, 
  sourceFormat, 
  targetFormat, 
  options,
  progressCallback,
  successCallback,
  errorCallback
) => {
  try {
    // Update progress
    progressCallback(20);
    
    // Create the conversion key
    const conversionKey = `${sourceFormat}_to_${targetFormat}`;
    
    // Map of conversion functions
    const conversions = {
      'fasta_to_genbank': fastaToGenbank,
      'fasta_to_embl': fastaToEmbl,
      'fasta_to_gff': fastaToGff,
      'fasta_to_clustal': fastaToClustal,
      'fasta_to_stockholm': fastaToStockholm,
      'fasta_to_pdb': fastaToPdb,
      'fastq_to_fasta': fastqToFasta,
      'genbank_to_fasta': genbankToFasta,
      'genbank_to_embl': genbankToEmbl,
      'genbank_to_gff': genbankToGff,
      'embl_to_fasta': emblToFasta,
      'embl_to_genbank': emblToGenbank,
      'embl_to_gff': emblToGff,
      'gff_to_fasta': gffToFasta,
      'gff_to_bed': gffToBed,
      'pdb_to_fasta': pdbToFasta,
      'pdb_to_mmcif': pdbToMmcif,
      'mmcif_to_pdb': mmcifToPdb,
      'mmcif_to_fasta': mmcifToFasta,
      'vcf_to_fasta': vcfToFasta,
      'newick_to_nexus': newickToNexus,
      'nexus_to_newick': nexusToNewick,
      'sam_to_fasta': samToFasta,
      'clustal_to_fasta': clustalToFasta,
      'clustal_to_stockholm': clustalToStockholm,
      'stockholm_to_fasta': stockholmToFasta,
      'stockholm_to_clustal': stockholmToClustal,
    };
    
    // Check if conversion is supported
    if (!conversions[conversionKey]) {
      throw new Error(`Conversion from ${sourceFormat} to ${targetFormat} is not supported`);
    }
    
    // Update progress
    progressCallback(40);
    
    // Perform the conversion
    const result = conversions[conversionKey](content, options);
    
    // Update progress
    progressCallback(80);
    
    // Return the result
    successCallback({
      content: result.content,
      format: targetFormat,
      originalFormat: sourceFormat,
      ...result
    });
    
  } catch (error) {
    errorCallback(`Conversion error: ${error.message}`);
  }
};

/**
 * Format a sequence with line breaks
 * @param {string} sequence - The sequence to format
 * @param {number} lineLength - The desired line length
 * @returns {string} The formatted sequence
 */
export const formatSequenceForFasta = (sequence, lineLength = 60) => {
  const lines = [];
  for (let i = 0; i < sequence.length; i += lineLength) {
    lines.push(sequence.substring(i, i + lineLength));
  }
  return lines.join('\n');
};

/**
 * Normalize a nucleotide sequence (uppercase and replace U with T)
 * @param {string} sequence - The sequence to normalize
 * @returns {string} The normalized sequence
 */
export const normalizeNucleotideSequence = (sequence) => {
  return sequence.toUpperCase().replace(/U/g, 'T');
};

/**
 * Remove gaps from a sequence
 * @param {string} sequence - The sequence with gaps
 * @returns {string} The sequence without gaps
 */
export const removeGapsFromSequence = (sequence) => {
  return sequence.replace(/-/g, '');
};

/**
 * Translate a nucleotide sequence to protein
 * @param {string} sequence - The nucleotide sequence
 * @returns {string} The protein sequence
 */
export const translateToProtein = (sequence) => {
  // Normalize the sequence first
  const normalizedSeq = normalizeNucleotideSequence(sequence);
  
  // Genetic code mapping
  const geneticCode = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
  };
  
  // Translate sequence
  let protein = '';
  for (let i = 0; i < normalizedSeq.length - 2; i += 3) {
    const codon = normalizedSeq.slice(i, i + 3);
    if (codon.length === 3) {
      protein += geneticCode[codon] || 'X'; // 'X' for unknown
    }
  }
  
  return protein;
};

/**
 * Extract sequences from FASTA content
 * @param {string} content - The FASTA content
 * @returns {Array} Array of sequence objects with headers and sequences
 */
export const extractFastaSequences = (content) => {
  const sequences = [];
  const lines = content.split('\n');
  
  let currentSeq = null;
  let currentHeader = "";
  
  for (const line of lines) {
    if (line.startsWith('>')) {
      if (currentSeq) {
        sequences.push({
          header: currentHeader,
          sequence: currentSeq
        });
      }
      currentHeader = line.substring(1).trim();
      currentSeq = "";
    } else if (line.trim() && currentHeader) {
      currentSeq += line.trim();
    }
  }
  
  if (currentSeq) {
    sequences.push({
      header: currentHeader,
      sequence: currentSeq
    });
  }
  
  return sequences;
};

/**
 * Convert FASTA to GenBank format
 * @param {string} content - The FASTA content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastaToGenbank = (content, options) => {
  const sequences = extractFastaSequences(content);
  let genbankContent = '';
  
  for (const seq of sequences) {
    // Generate a locus name from the header
    const locusParts = seq.header.split(/\s+/);
    const locusName = locusParts[0].replace(/[^a-zA-Z0-9]/g, '_').substring(0, 15);
    
    // Get sequence length
    const seqLength = seq.sequence.length;
    
    // Format the sequence for GenBank (10 bases per group, 6 groups per line)
    const formattedSeq = [];
    for (let i = 0; i < seq.sequence.length; i += 60) {
      let line = '';
      let lineSeq = seq.sequence.substring(i, i + 60);
      
      // Right-align the base position number
      let basePos = i + 1;
      line += ' '.repeat(9 - basePos.toString().length) + basePos;
      
      // Format the sequence in groups of 10
      for (let j = 0; j < lineSeq.length; j += 10) {
        line += ' ' + lineSeq.substring(j, j + 10);
      }
      
      formattedSeq.push(line);
    }
    
    // Create the GenBank record
    genbankContent += `LOCUS       ${locusName.padEnd(15)} ${seqLength} bp    DNA     linear   01-JAN-2023
DEFINITION  ${seq.header}
ACCESSION   
VERSION     
SOURCE      .
  ORGANISM  .
FEATURES             Location/Qualifiers
     source          1..${seqLength}
ORIGIN
${formattedSeq.join('\n')}
//

`;
  }
  
  return {
    content: genbankContent,
    sequenceCount: sequences.length,
    totalLength: sequences.reduce((sum, seq) => sum + seq.sequence.length, 0),
    note: 'Basic conversion without annotations. Locus names were derived from sequence headers.'
  };
};

/**
 * Convert FASTA to EMBL format
 * @param {string} content - The FASTA content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastaToEmbl = (content, options) => {
  const sequences = extractFastaSequences(content);
  let emblContent = '';
  
  for (const seq of sequences) {
    // Generate ID from header
    const idParts = seq.header.split(/\s+/);
    const idName = idParts[0].replace(/[^a-zA-Z0-9]/g, '_').substring(0, 15);
    
    // Get sequence length
    const seqLength = seq.sequence.length;
    
    // Format the sequence for EMBL (10 bases per group, 6 groups per line)
    const formattedSeq = [];
    for (let i = 0; i < seq.sequence.length; i += 60) {
      let lineSeq = seq.sequence.substring(i, i + 60);
      let groups = [];
      
      for (let j = 0; j < lineSeq.length; j += 10) {
        groups.push(lineSeq.substring(j, j + 10));
      }
      
      let line = '     ' + groups.join(' ') + ' ' + (i + 60);
      formattedSeq.push(line);
    }
    
    // Create the EMBL record
    emblContent += `ID   ${idName}; SV 1; linear; DNA; STD; UNC; ${seqLength} BP.
XX
AC   ;
XX
DE   ${seq.header}
XX
FH   Key             Location/Qualifiers
FT   source          1..${seqLength}
XX
SQ   Sequence ${seqLength} BP;
${formattedSeq.join('\n')}
//

`;
  }
  
  return {
    content: emblContent,
    sequenceCount: sequences.length,
    totalLength: sequences.reduce((sum, seq) => sum + seq.sequence.length, 0),
    note: 'Basic conversion without annotations. IDs were derived from sequence headers.'
  };
};

/**
 * Convert FASTA to GFF format
 * @param {string} content - The FASTA content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastaToGff = (content, options) => {
  const sequences = extractFastaSequences(content);
  let gffContent = '##gff-version 3\n';
  
  for (const seq of sequences) {
    // Get sequence ID from header
    const seqId = seq.header.split(/\s+/)[0];
    
    // Get sequence length
    const seqLength = seq.sequence.length;
    
    // Add sequence-region pragma
    gffContent += `##sequence-region ${seqId} 1 ${seqLength}\n`;
    
    // Add a basic 'region' feature for each sequence
    gffContent += `${seqId}\t.\tregion\t1\t${seqLength}\t.\t+\t.\tID=${seqId};Name=${seqId}\n`;
    
    // If requested, we could add CDS, exon, or other features here
  }
  
  return {
    content: gffContent,
    sequenceCount: sequences.length,
    totalLength: sequences.reduce((sum, seq) => sum + seq.sequence.length, 0),
    note: 'Created minimal GFF with region features. No gene or CDS predictions were performed.'
  };
};

/**
 * Convert FASTA to Clustal format
 * @param {string} content - The FASTA content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastaToClustal = (content, options) => {
  const sequences = extractFastaSequences(content);
  
  if (sequences.length < 2) {
    throw new Error('CLUSTAL format requires at least 2 sequences');
  }
  
  // Extract sequence IDs and prepare data
  const seqData = sequences.map(seq => {
    return {
      id: seq.header.split(/\s+/)[0],
      sequence: seq.sequence
    };
  });
  
  // Find max sequence ID length for padding
  const maxIdLength = Math.max(...seqData.map(seq => seq.id.length));
  
  // Create CLUSTAL header
  let clustalContent = 'CLUSTAL W (1.82) multiple sequence alignment\n\n';
  
  // Format sequences in blocks of 60 characters
  const lineLength = 60;
  const maxLength = Math.max(...seqData.map(seq => seq.sequence.length));
  
  for (let i = 0; i < maxLength; i += lineLength) {
    // Add a blank line between blocks
    if (i > 0) {
      clustalContent += '\n';
    }
    
    // Add each sequence line
    for (const seq of seqData) {
      const segment = seq.sequence.substring(i, i + lineLength);
      if (segment.length > 0) {
        clustalContent += `${seq.id.padEnd(maxIdLength)} ${segment}\n`;
      }
    }
    
    // Add conservation line (simplistic)
    clustalContent += ' '.repeat(maxIdLength) + ' ';
    const firstSeq = seqData[0].sequence.substring(i, i + lineLength);
    
    for (let j = 0; j < firstSeq.length; j++) {
      // Check if all sequences have the same character at this position
      const allMatch = seqData.every(seq => 
        j < seq.sequence.length - i && seq.sequence[i + j] === firstSeq[j]
      );
      
      clustalContent += allMatch ? '*' : ' ';
    }
    
    clustalContent += '\n';
  }
  
  return {
    content: clustalContent,
    sequenceCount: sequences.length,
    totalLength: sequences.reduce((sum, seq) => sum + seq.sequence.length, 0),
    note: 'Converted to CLUSTAL without performing alignment. If these sequences are not pre-aligned, the result may not be meaningful.'
  };
};

/**
 * Convert FASTA to Stockholm format
 * @param {string} content - The FASTA content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastaToStockholm = (content, options) => {
  const sequences = extractFastaSequences(content);
  
  if (sequences.length < 2) {
    throw new Error('Stockholm format is designed for multiple sequence alignments');
  }
  
  // Extract sequence IDs and prepare data
  const seqData = sequences.map(seq => {
    return {
      id: seq.header.split(/\s+/)[0],
      sequence: seq.sequence
    };
  });
  
  // Find max sequence ID length for padding
  const maxIdLength = Math.max(...seqData.map(seq => seq.id.length));
  
  // Create Stockholm header
  let stockholmContent = '# STOCKHOLM 1.0\n';
  
  // Add sequences
  for (const seq of seqData) {
    stockholmContent += `${seq.id.padEnd(maxIdLength)} ${seq.sequence}\n`;
  }
  
  // Add a simple consensus line
  stockholmContent += '#=GC RF' + ' '.repeat(maxIdLength - 6) + ' ';
  
  // Generate a simple consensus (using * for identical positions)
  const firstSeq = seqData[0].sequence;
  for (let i = 0; i < firstSeq.length; i++) {
    const allMatch = seqData.every(seq => 
      i < seq.sequence.length && seq.sequence[i] === firstSeq[i]
    );
    
    stockholmContent += allMatch ? '*' : '.';
  }
  
  // Close the Stockholm file
  stockholmContent += '\n//\n';
  
  return {
    content: stockholmContent,
    sequenceCount: sequences.length,
    totalLength: sequences.reduce((sum, seq) => sum + seq.sequence.length, 0),
    note: 'Converted to Stockholm without performing alignment. If these sequences are not pre-aligned, the result may not be meaningful.'
  };
};

/**
 * Convert FASTA to PDB format (basic conversion for protein sequences)
 * @param {string} content - The FASTA content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastaToPdb = (content, options) => {
  const sequences = extractFastaSequences(content);
  
  if (sequences.length > 1) {
    throw new Error('This converter can only convert a single protein sequence to a basic PDB file');
  }
  
  const sequence = sequences[0];
  const seqId = sequence.header.split(/\s+/)[0];
  
  // Basic PDB format with just CA atoms in a linear chain
  let pdbContent = `HEADER    PROTEIN                                 01-JAN-23   ${seqId}\n`;
  pdbContent += `TITLE     CONVERTED FROM FASTA\n`;
  pdbContent += `REMARK    This is a simplified PDB file generated from a FASTA sequence\n`;
  pdbContent += `REMARK    It contains only alpha carbon atoms in a linear arrangement\n`;
  
  // Create CA atoms for each residue (pseudo-structure)
  let atomNum = 1;
  const chainId = 'A';
  
  for (let i = 0; i < sequence.sequence.length; i++) {
    const residue = sequence.sequence[i];
    const residueNum = i + 1;
    
    // Simple linear placement of CA atoms along the x-axis
    const x = i * 3.8; // ~3.8 Å is typical CA-CA distance
    const y = 0.0;
    const z = 0.0;
    
    // Format PDB ATOM record
    pdbContent += `ATOM  ${atomNum.toString().padStart(5)} ${' CA '.padEnd(4)} ${residue} ${chainId}${residueNum.toString().padStart(4)}    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00           C\n`;
    
    atomNum++;
  }
  
  pdbContent += `TER   ${atomNum.toString().padStart(5)}      ${sequence.sequence[sequence.sequence.length - 1]} ${chainId}${sequence.sequence.length.toString().padStart(4)}\n`;
  pdbContent += 'END\n';
  
  return {
    content: pdbContent,
    sequenceCount: 1,
    totalLength: sequence.sequence.length,
    note: 'Created a minimal PDB file with only CA atoms in a linear arrangement. This is not a real 3D structure.'
  };
};

/**
 * Convert FASTQ to FASTA format
 * @param {string} content - The FASTQ content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const fastqToFasta = (content, options) => {
  const lines = content.split('\n');
  let fastaContent = '';
  let sequenceCount = 0;
  let totalLength = 0;
  
  // Process FASTQ in groups of 4 lines
  for (let i = 0; i < lines.length; i += 4) {
    if (i + 3 < lines.length && lines[i].startsWith('@')) {
      const header = lines[i].substring(1); // Remove @ and use the header
      const sequence = lines[i + 1];
      
      if (sequence) {
        // Add to FASTA format
        if (options.formatSequence) {
          fastaContent += `>${header}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
        } else {
          fastaContent += `>${header}\n${sequence}\n`;
        }
        
        sequenceCount++;
        totalLength += sequence.length;
      }
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount,
    totalLength,
    note: options.includeQuality ? 'Quality scores were discarded during conversion to FASTA.' : undefined
  };
};

/**
 * Convert GenBank to FASTA format
 * @param {string} content - The GenBank content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const genbankToFasta = (content, options) => {
  const sections = content.split('//');
  let fastaContent = '';
  let sequenceCount = 0;
  let totalLength = 0;
  
  for (const section of sections) {
    if (!section.trim()) continue;
    
    // Extract LOCUS information
    const locusMatch = section.match(/LOCUS\s+(\S+)/);
    const locus = locusMatch ? locusMatch[1] : 'Unknown';
    
    // Extract DEFINITION
    const defMatch = section.match(/DEFINITION\s+(.*?)(?=\n\S)/s);
    const definition = defMatch ? defMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract ACCESSION
    const accMatch = section.match(/ACCESSION\s+(\S+)/);
    const accession = accMatch ? accMatch[1] : '';
    
    // Extract ORGANISM
    const orgMatch = section.match(/\s+ORGANISM\s+(.*?)(?=\n\S)/s);
    const organism = orgMatch ? orgMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract sequence
    const originMatch = section.match(/ORIGIN([\s\S]*?)(?=$)/);
    if (originMatch) {
      let sequence = originMatch[1].replace(/\d+|\s+/g, '');
      
      if (sequence) {
        // Create a descriptive FASTA header
        let header = `${locus}`;
        if (accession) header += ` ${accession}`;
        if (definition) header += ` ${definition}`;
        if (organism) header += ` [${organism}]`;
        
        if (options.formatSequence) {
          fastaContent += `>${header}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
        } else {
          fastaContent += `>${header}\n${sequence}\n`;
        }
        
        sequenceCount++;
        totalLength += sequence.length;
      }
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount,
    totalLength,
    note: 'Features and annotations were not included in the FASTA conversion.'
  };
};

/**
 * Convert GenBank to EMBL format
 * @param {string} content - The GenBank content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const genbankToEmbl = (content, options) => {
  const sections = content.split('//');
  let emblContent = '';
  let sequenceCount = 0;
  let totalLength = 0;
  
  for (const section of sections) {
    if (!section.trim()) continue;
    
    // Extract LOCUS information
    const locusMatch = section.match(/LOCUS\s+(\S+)\s+(\d+)\s+bp/);
    const locus = locusMatch ? locusMatch[1] : 'Unknown';
    const length = locusMatch ? parseInt(locusMatch[2]) : 0;
    
    // Extract DEFINITION
    const defMatch = section.match(/DEFINITION\s+(.*?)(?=\n\S)/s);
    const definition = defMatch ? defMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract ACCESSION
    const accMatch = section.match(/ACCESSION\s+(\S+)/);
    const accession = accMatch ? accMatch[1] : '';
    
    // Extract ORGANISM
    const orgMatch = section.match(/\s+ORGANISM\s+(.*?)(?=\n\S)/s);
    const organism = orgMatch ? orgMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract FEATURES
    const featuresMatch = section.match(/FEATURES\s+Location\/Qualifiers\s+(.*?)(?=ORIGIN|\n\S\S)/s);
    let features = '';
    
    if (featuresMatch) {
      const featureBlock = featuresMatch[1];
      const featureLines = featureBlock.split('\n');
      
      // Convert GenBank features to EMBL format
      for (const line of featureLines) {
        if (line.match(/^\s{5}\S/)) {
          // Feature line
          const parts = line.trim().split(/\s+/, 2);
          features += `FT   ${parts[0].padEnd(15)}${parts[1] || ''}\n`;
        } else if (line.match(/^\s{21}\/\S/)) {
          // Qualifier line, convert /qualifier="value" to /qualifier=value
          features += line.replace(/^\s{21}/, 'FT                   ') + '\n';
        } else if (line.trim()) {
          // Continuation of a qualifier
          features += line.replace(/^\s{21}/, 'FT                   ') + '\n';
        }
      }
    }
    
    // Extract sequence
    const originMatch = section.match(/ORIGIN([\s\S]*?)(?=$)/);
    let sequence = '';
    
    if (originMatch) {
      sequence = originMatch[1].replace(/\d+|\s+/g, '');
      totalLength += sequence.length;
      sequenceCount++;
      
      // Format the sequence for EMBL (10 bases per group, 6 groups per line)
      const formattedSeq = [];
      for (let i = 0; i < sequence.length; i += 60) {
        let lineSeq = sequence.substring(i, i + 60);
        let groups = [];
        
        for (let j = 0; j < lineSeq.length; j += 10) {
          groups.push(lineSeq.substring(j, j + 10));
        }
        
        let line = '     ' + groups.join(' ') + ' ' + (i + 60);
        formattedSeq.push(line);
      }
      
      // Create the EMBL record
      emblContent += `ID   ${locus}; SV 1; linear; DNA; STD; UNC; ${length} BP.\n`;
      emblContent += `XX\n`;
      
      if (accession) {
        emblContent += `AC   ${accession};\n`;
        emblContent += `XX\n`;
      }
      
      if (definition) {
        emblContent += `DE   ${definition}\n`;
        emblContent += `XX\n`;
      }
      
      if (organism) {
        emblContent += `OS   ${organism}\n`;
        emblContent += `XX\n`;
      }
      
      if (features) {
        emblContent += `FH   Key             Location/Qualifiers\n`;
        emblContent += features;
        emblContent += `XX\n`;
      }
      
      emblContent += `SQ   Sequence ${length} BP;\n`;
      emblContent += formattedSeq.join('\n') + '\n';
      emblContent += `//\n\n`;
    }
  }
  
  return {
    content: emblContent,
    sequenceCount,
    totalLength,
    note: 'Some GenBank-specific annotations might not have direct EMBL equivalents.'
  };
};

/**
 * Convert GenBank to GFF format
 * @param {string} content - The GenBank content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const genbankToGff = (content, options) => {
  const sections = content.split('//');
  let gffContent = '##gff-version 3\n';
  let sequenceCount = 0;
  let featureCount = 0;
  
  for (const section of sections) {
    if (!section.trim()) continue;
    
    // Extract LOCUS information
    const locusMatch = section.match(/LOCUS\s+(\S+)\s+(\d+)\s+bp/);
    const locus = locusMatch ? locusMatch[1] : 'Unknown';
    const length = locusMatch ? parseInt(locusMatch[2]) : 0;
    
    // Add sequence-region pragma
    gffContent += `##sequence-region ${locus} 1 ${length}\n`;
    
    // Extract FEATURES
    const featuresMatch = section.match(/FEATURES\s+Location\/Qualifiers\s+(.*?)(?=ORIGIN|\n\S\S)/s);
    
    if (featuresMatch) {
      const featureBlock = featuresMatch[1];
      const featureLines = featureBlock.split('\n');
      
      let currentFeature = null;
      let currentAttributes = '';
      
      for (let i = 0; i < featureLines.length; i++) {
        const line = featureLines[i];
        
        if (line.match(/^\s{5}\S/)) {
          // Save previous feature
          if (currentFeature) {
            gffContent += `${locus}\tGenBank\t${currentFeature.type}\t${currentFeature.start}\t${currentFeature.end}\t.\t${currentFeature.strand}\t.\t${currentAttributes}\n`;
            featureCount++;
          }
          
          // New feature
          const parts = line.trim().split(/\s+/, 2);
          const type = parts[0];
          const location = parts[1] || '';
          
          // Parse location
          let start = 1;
          let end = length;
          let strand = '+';
          
          // Parse common location formats
          const simpleRangeMatch = location.match(/(\d+)\.\.(\d+)/);
          if (simpleRangeMatch) {
            start = parseInt(simpleRangeMatch[1]);
            end = parseInt(simpleRangeMatch[2]);
          }
          
          const revComplementMatch = location.match(/complement\((\d+)\.\.(\d+)\)/);
          if (revComplementMatch) {
            start = parseInt(revComplementMatch[1]);
            end = parseInt(revComplementMatch[2]);
            strand = '-';
          }
          
          // Handle other types of locations
          if (location.match(/^\d+$/)) {
            // Single base feature
            start = parseInt(location);
            end = start;
          }
          
          if (location.includes('join')) {
            // For joined features, use the start of the first and end of the last
            const joinMatch = location.match(/join\((.*)\)/);
            if (joinMatch) {
              const parts = joinMatch[1].split(',');
              const firstPart = parts[0].match(/\d+/);
              const lastPart = parts[parts.length - 1].match(/\d+\.\.(\d+)/);
              
              if (firstPart && lastPart) {
                start = parseInt(firstPart[0]);
                end = parseInt(lastPart[1]);
              }
            }
          }
          
          currentFeature = { type, start, end, strand };
          currentAttributes = `ID=${locus}_${type}_${featureCount + 1};`;
        } else if (line.match(/^\s{21}\/\S/) && currentFeature) {
          // Qualifier line
          const match = line.trim().match(/\/(\S+)=(?:"([^"]*)"|([\S]+))/);
          if (match) {
            const key = match[1];
            const value = match[2] || match[3];
            
            // Convert GenBank qualifiers to GFF attributes
            if (key === 'gene') {
              currentAttributes += `Name=${value};`;
            } else if (key === 'locus_tag') {
              currentAttributes += `locus_tag=${value};`;
            } else if (key === 'product') {
              currentAttributes += `product=${value};`;
            } else if (key === 'note') {
              currentAttributes += `Note=${value};`;
            } else if (key === 'db_xref') {
              currentAttributes += `Dbxref=${value};`;
            } else {
              currentAttributes += `${key}=${value};`;
            }
          }
        }
      }
      
      // Don't forget the last feature
      if (currentFeature) {
        gffContent += `${locus}\tGenBank\t${currentFeature.type}\t${currentFeature.start}\t${currentFeature.end}\t.\t${currentFeature.strand}\t.\t${currentAttributes}\n`;
        featureCount++;
      }
    }
    
    sequenceCount++;
  }
  
  return {
    content: gffContent,
    sequenceCount,
    featureCount,
    note: 'Some complex GenBank locations may not convert perfectly to GFF.'
  };
};

/**
 * Convert EMBL to FASTA format
 * @param {string} content - The EMBL content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const emblToFasta = (content, options) => {
  const sections = content.split('//');
  let fastaContent = '';
  let sequenceCount = 0;
  let totalLength = 0;
  
  for (const section of sections) {
    if (!section.trim()) continue;
    
    // Extract ID information
    const idMatch = section.match(/ID\s+([^;]+)/);
    const id = idMatch ? idMatch[1].trim() : 'Unknown';
    
    // Extract AC (accession)
    const acMatch = section.match(/AC\s+([^;]+)/);
    const accession = acMatch ? acMatch[1].trim() : '';
    
    // Extract DE (description)
    const deMatch = section.match(/DE\s+(.*?)(?=\nXX|\n\S\S)/s);
    const description = deMatch ? deMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract OS (organism)
    const osMatch = section.match(/OS\s+(.*?)(?=\nXX|\n\S\S)/s);
    const organism = osMatch ? osMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract sequence
    const sqMatch = section.match(/SQ.*?\n([\s\S]*?)(?=\/\/|$)/s);
    
    if (sqMatch) {
      let sequence = sqMatch[1].replace(/\d+|\s+/g, '');
      
      if (sequence) {
        // Create a descriptive FASTA header
        let header = `${id}`;
        if (accession) header += ` ${accession}`;
        if (description) header += ` ${description}`;
        if (organism) header += ` [${organism}]`;
        
        if (options.formatSequence) {
          fastaContent += `>${header}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
        } else {
          fastaContent += `>${header}\n${sequence}\n`;
        }
        
        sequenceCount++;
        totalLength += sequence.length;
      }
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount,
    totalLength,
    note: 'Features and annotations were not included in the FASTA conversion.'
  };
};

/**
 * Convert EMBL to GenBank format
 * @param {string} content - The EMBL content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const emblToGenbank = (content, options) => {
  const sections = content.split('//');
  let genbankContent = '';
  let sequenceCount = 0;
  let totalLength = 0;
  
  for (const section of sections) {
    if (!section.trim()) continue;
    
    // Extract ID information
    const idMatch = section.match(/ID\s+([^;]+);\s+SV\s+(\d+);\s+.*?;\s+(\d+)\s+BP/);
    const id = idMatch ? idMatch[1].trim() : 'Unknown';
    const version = idMatch ? idMatch[2] : '1';
    const length = idMatch ? parseInt(idMatch[3]) : 0;
    
    // Extract AC (accession)
    const acMatch = section.match(/AC\s+([^;]+)/);
    const accession = acMatch ? acMatch[1].trim() : '';
    
    // Extract DE (description)
    const deMatch = section.match(/DE\s+(.*?)(?=\nXX|\n\S\S)/s);
    const description = deMatch ? deMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract OS (organism)
    const osMatch = section.match(/OS\s+(.*?)(?=\nXX|\n\S\S)/s);
    const organism = osMatch ? osMatch[1].replace(/\n\s+/g, ' ').trim() : '';
    
    // Extract features
    const featureBlock = section.match(/FH.*?\nFT([\s\S]*?)(?=XX|$)/s);
    let features = '';
    
    if (featureBlock) {
      const featureLines = featureBlock[1].split('\n');
      
      // Convert EMBL features to GenBank format
      for (const line of featureLines) {
        if (line.match(/^\s+\S+\s+/)) {
          // Feature line
          const match = line.match(/^\s+(\S+)\s+(.*)/);
          if (match) {
            const type = match[1];
            const location = match[2];
            features += `     ${type.padEnd(15)}${location}\n`;
          }
        } else if (line.match(/^\s+\/\S/)) {
          // Qualifier line
          features += line.replace(/^\s+\//, '                     /') + '\n';
        }
      }
    }
    
    // Extract sequence
    const sqMatch = section.match(/SQ.*?\n([\s\S]*?)(?=\/\/|$)/s);
    let sequence = '';
    
    if (sqMatch) {
      sequence = sqMatch[1].replace(/\d+|\s+/g, '');
      totalLength += sequence.length;
      sequenceCount++;
      
      // Format the sequence for GenBank (10 bases per group, 6 groups per line)
      const formattedSeq = [];
      for (let i = 0; i < sequence.length; i += 60) {
        let line = '';
        let lineSeq = sequence.substring(i, i + 60);
        
        // Right-align the base position number
        let basePos = i + 1;
        line += ' '.repeat(9 - basePos.toString().length) + basePos;
        
        // Format the sequence in groups of 10
        for (let j = 0; j < lineSeq.length; j += 10) {
          line += ' ' + lineSeq.substring(j, j + 10);
        }
        
        formattedSeq.push(line);
      }
      
      // Create the GenBank record
      genbankContent += `LOCUS       ${id.padEnd(15)} ${length} bp    DNA     linear   01-JAN-2023\n`;
      
      if (description) {
        genbankContent += `DEFINITION  ${description}\n`;
      }
      
      if (accession) {
        genbankContent += `ACCESSION   ${accession}\n`;
        genbankContent += `VERSION     ${accession}.${version}\n`;
      }
      
      genbankContent += `SOURCE      .\n`;
      
      if (organism) {
        genbankContent += `  ORGANISM  ${organism}\n`;
        genbankContent += `            .\n`;
      }
      
      if (features) {
        genbankContent += `FEATURES             Location/Qualifiers\n`;
        genbankContent += features;
      }
      
      genbankContent += `ORIGIN\n`;
      genbankContent += formattedSeq.join('\n') + '\n';
      genbankContent += `//\n\n`;
    }
  }
  
  return {
    content: genbankContent,
    sequenceCount,
    totalLength,
    note: 'Some EMBL-specific annotations might not have direct GenBank equivalents.'
  };
};

/**
 * Convert EMBL to GFF format
 * @param {string} content - The EMBL content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const emblToGff = (content, options) => {
  const sections = content.split('//');
  let gffContent = '##gff-version 3\n';
  let sequenceCount = 0;
  let featureCount = 0;
  
  for (const section of sections) {
    if (!section.trim()) continue;
    
    // Extract ID information
    const idMatch = section.match(/ID\s+([^;]+);\s+SV\s+(\d+);\s+.*?;\s+(\d+)\s+BP/);
    const id = idMatch ? idMatch[1].trim() : 'Unknown';
    const length = idMatch ? parseInt(idMatch[3]) : 0;
    
    // Add sequence-region pragma
    gffContent += `##sequence-region ${id} 1 ${length}\n`;
    
    // Extract features
    const featureBlock = section.match(/FH.*?\nFT([\s\S]*?)(?=XX|$)/s);
    
    if (featureBlock) {
      const featureLines = featureBlock[1].split('\n');
      
      let currentFeature = null;
      let currentAttributes = '';
      
      for (let i = 0; i < featureLines.length; i++) {
        const line = featureLines[i];
        
        // Feature line
        if (line.match(/^\s+\S+\s+\S/)) {
          // Save previous feature
          if (currentFeature) {
            gffContent += `${id}\tEMBL\t${currentFeature.type}\t${currentFeature.start}\t${currentFeature.end}\t.\t${currentFeature.strand}\t.\t${currentAttributes}\n`;
            featureCount++;
          }
          
          // New feature
          const match = line.match(/^\s+(\S+)\s+(.*)/);
          if (match) {
            const type = match[1];
            const location = match[2];
            
            // Parse location
            let start = 1;
            let end = length;
            let strand = '+';
            
            // Parse common location formats
            const simpleRangeMatch = location.match(/(\d+)\.\.(\d+)/);
            if (simpleRangeMatch) {
              start = parseInt(simpleRangeMatch[1]);
              end = parseInt(simpleRangeMatch[2]);
            }
            
            const revComplementMatch = location.match(/complement\((\d+)\.\.(\d+)\)/);
            if (revComplementMatch) {
              start = parseInt(revComplementMatch[1]);
              end = parseInt(revComplementMatch[2]);
              strand = '-';
            }
            
            // Handle other types of locations
            if (location.match(/^\d+$/)) {
              // Single base feature
              start = parseInt(location);
              end = start;
            }
            
            if (location.includes('join')) {
              // For joined features, use the start of the first and end of the last
              const joinMatch = location.match(/join\((.*)\)/);
              if (joinMatch) {
                const parts = joinMatch[1].split(',');
                const firstPart = parts[0].match(/\d+/);
                const lastPart = parts[parts.length - 1].match(/\d+\.\.(\d+)/);
                
                if (firstPart && lastPart) {
                  start = parseInt(firstPart[0]);
                  end = parseInt(lastPart[1]);
                }
              }
            }
            
            currentFeature = { type, start, end, strand };
            currentAttributes = `ID=${id}_${type}_${featureCount + 1};`;
          }
        } else if (line.match(/^\s+\/\S/) && currentFeature) {
          // Qualifier line
          const match = line.trim().match(/\/(\S+)=(?:"([^"]*)"|([\S]+))/);
          if (match) {
            const key = match[1];
            const value = match[2] || match[3];
            
            // Convert EMBL qualifiers to GFF attributes
            if (key === 'gene') {
              currentAttributes += `Name=${value};`;
            } else if (key === 'locus_tag') {
              currentAttributes += `locus_tag=${value};`;
            } else if (key === 'product') {
              currentAttributes += `product=${value};`;
            } else if (key === 'note') {
              currentAttributes += `Note=${value};`;
            } else if (key === 'db_xref') {
              currentAttributes += `Dbxref=${value};`;
            } else {
              currentAttributes += `${key}=${value};`;
            }
          }
        }
      }
      
      // Don't forget the last feature
      if (currentFeature) {
        gffContent += `${id}\tEMBL\t${currentFeature.type}\t${currentFeature.start}\t${currentFeature.end}\t.\t${currentFeature.strand}\t.\t${currentAttributes}\n`;
        featureCount++;
      }
    }
    
    sequenceCount++;
  }
  
  return {
    content: gffContent,
    sequenceCount,
    featureCount,
    note: 'Some complex EMBL locations may not convert perfectly to GFF.'
  };
};

/**
 * Convert GFF to FASTA format (mock implementation without actual sequence data)
 * @param {string} content - The GFF content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const gffToFasta = (content, options) => {
  // This is a placeholder as GFF files don't typically contain sequence data
  // In a real implementation, this would require the reference sequence as additional input
  
  const lines = content.split('\n');
  const sequenceRegions = [];
  
  // Extract sequence-region pragmas
  for (const line of lines) {
    if (line.startsWith('##sequence-region')) {
      const parts = line.split(/\s+/);
      if (parts.length >= 4) {
        sequenceRegions.push({
          id: parts[1],
          start: parseInt(parts[2]),
          end: parseInt(parts[3]),
          length: parseInt(parts[3]) - parseInt(parts[2]) + 1
        });
      }
    }
  }
  
  let fastaContent = '';
  let totalLength = 0;
  
  // Create placeholder sequences based on sequence-region pragmas
  for (const region of sequenceRegions) {
    fastaContent += `>${region.id} (placeholder sequence for GFF region)\n`;
    
    // Create a dummy sequence of Ns
    const dummySequence = 'N'.repeat(region.length);
    
    if (options.formatSequence) {
      fastaContent += formatSequenceForFasta(dummySequence, options.lineLength || 60) + '\n';
    } else {
      fastaContent += dummySequence + '\n';
    }
    
    totalLength += region.length;
  }
  
  if (sequenceRegions.length === 0) {
    fastaContent = '>GFF_derived_sequence (placeholder)\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n';
  }
  
  return {
    content: fastaContent,
    sequenceCount: sequenceRegions.length || 1,
    totalLength: totalLength || 60,
    note: 'GFF files typically do not contain sequence data. This FASTA contains placeholder sequences only. For a complete conversion, provide the reference sequence file separately.'
  };
};

/**
 * Convert GFF to BED format
 * @param {string} content - The GFF content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const gffToBed = (content, options) => {
  const lines = content.split('\n');
  let bedContent = '';
  let featureCount = 0;
  
  for (const line of lines) {
    // Skip comment and empty lines
    if (line.startsWith('#') || !line.trim()) {
      continue;
    }
    
    // Parse GFF fields
    const fields = line.split('\t');
    if (fields.length < 8) {
      continue;
    }
    
    const seqId = fields[0];
    const source = fields[1];
    const type = fields[2];
    const start = parseInt(fields[3]);
    const end = parseInt(fields[4]);
    const score = fields[5] === '.' ? '0' : fields[5];
    const strand = fields[6] === '.' ? '+' : fields[6];
    const phase = fields[7];
    const attributes = fields[8] || '';
    
    // Parse attributes to extract name and other fields
    let name = `${type}_${featureCount + 1}`;
    const idMatch = attributes.match(/ID=([^;]+)/);
    const nameMatch = attributes.match(/Name=([^;]+)/);
    
    if (nameMatch) {
      name = nameMatch[1];
    } else if (idMatch) {
      name = idMatch[1];
    }
    
    // Convert to BED (remember BED is 0-based for start, 1-based for end)
    bedContent += `${seqId}\t${start - 1}\t${end}\t${name}\t${score}\t${strand}\n`;
    
    featureCount++;
  }
  
  return {
    content: bedContent,
    featureCount,
    note: 'GFF attributes were simplified to fit BED format. Only ID or Name attributes were used for the name field.'
  };
};

/**
 * Convert PDB to FASTA format (extract protein sequence)
 * @param {string} content - The PDB content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const pdbToFasta = (content, options) => {
  const lines = content.split('\n');
  const chains = new Map();
  
  // Amino acid 3-letter to 1-letter conversion
  const aa3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J',
    'XAA': 'X', 'UNK': 'X'
  };
  
  // Extract header info for metadata
  let headerInfo = '';
  const headerMatch = content.match(/HEADER\s+(.*)/);
  const titleMatch = content.match(/TITLE\s+(.*?)(?=\nREMARK|\nSOURCE|\nAUTHOR|\nATOM)/s);
  
  if (headerMatch) {
    headerInfo += headerMatch[1].trim() + ' ';
  }
  
  if (titleMatch) {
    headerInfo += titleMatch[1].replace(/\n\s+/g, ' ').trim();
  }
  
  // Extract protein sequences from ATOM records
  for (const line of lines) {
    if (line.startsWith('ATOM')) {
      const residueName = line.substring(17, 20).trim();
      const chainId = line.substring(21, 22).trim();
      const residueNum = parseInt(line.substring(22, 26).trim());
      const atomName = line.substring(12, 16).trim();
      
      // Only process CA atoms (alpha carbon) to avoid duplicates
      if (atomName === 'CA' && residueName in aa3to1) {
        if (!chains.has(chainId)) {
          chains.set(chainId, new Map());
        }
        
        const chain = chains.get(chainId);
        chain.set(residueNum, aa3to1[residueName]);
      }
    }
  }
  
  // Convert to FASTA format
  let fastaContent = '';
  let totalLength = 0;
  
  for (const [chainId, residues] of chains.entries()) {
    // Sort residues by residue number
    const sortedResidues = new Map([...residues.entries()].sort((a, b) => a[0] - b[0]));
    
    // Construct sequence
    const sequence = [...sortedResidues.values()].join('');
    totalLength += sequence.length;
    
    // Create header
    const header = headerInfo ? 
      `${headerInfo} Chain ${chainId}` : 
      `PDB Chain ${chainId}`;
    
    if (options.formatSequence) {
      fastaContent += `>${header}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
    } else {
      fastaContent += `>${header}\n${sequence}\n`;
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount: chains.size,
    totalLength,
    note: 'Sequences were extracted from ATOM records for standard amino acids only.'
  };
};

/**
 * Convert PDB to mmCIF format (simplified implementation)
 * @param {string} content - The PDB content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const pdbToMmcif = (content, options) => {
  const lines = content.split('\n');
  let mmcifContent = 'data_structure\n';
  
  // Extract header info
  let title = 'Unknown';
  const titleMatch = content.match(/TITLE\s+(.*?)(?=\nREMARK|\nSOURCE|\nAUTHOR|\nATOM)/s);
  if (titleMatch) {
    title = titleMatch[1].replace(/\n\s+/g, ' ').trim();
  }
  
  // Add structure metadata
  mmcifContent += `_entry.id structure\n`;
  mmcifContent += `_struct.title "${title}"\n\n`;
  
  // Collect atom data
  const atoms = [];
  const models = new Set();
  const chains = new Set();
  const entities = new Map();
  
  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const recordType = line.substring(0, 6).trim();
      const atomNum = parseInt(line.substring(6, 11).trim());
      const atomName = line.substring(12, 16).trim();
      const altLoc = line.substring(16, 17).trim();
      const residueName = line.substring(17, 20).trim();
      const chainId = line.substring(21, 22).trim();
      const residueNum = parseInt(line.substring(22, 26).trim());
      const iCode = line.substring(26, 27).trim();
      const x = parseFloat(line.substring(30, 38).trim());
      const y = parseFloat(line.substring(38, 46).trim());
      const z = parseFloat(line.substring(46, 54).trim());
      const occupancy = parseFloat(line.substring(54, 60).trim() || "1.0");
      const tempFactor = parseFloat(line.substring(60, 66).trim() || "0.0");
      const element = line.substring(76, 78).trim();
      
      // Add to collections
      atoms.push({
        recordType, atomNum, atomName, altLoc, residueName, chainId, 
        residueNum, iCode, x, y, z, occupancy, tempFactor, element
      });
      
      chains.add(chainId);
      
      if (!entities.has(chainId)) {
        entities.set(chainId, new Set());
      }
      entities.get(chainId).add(residueName);
    } else if (line.startsWith('MODEL')) {
      const modelNum = parseInt(line.substring(10, 14).trim());
      models.add(modelNum);
    }
  }
  
  // Define entities
  mmcifContent += '# Entity information\n';
  mmcifContent += 'loop_\n';
  mmcifContent += '_entity.id\n';
  mmcifContent += '_entity.type\n';
  mmcifContent += '_entity.pdbx_description\n';
  
  let entityId = 1;
  for (const [chainId, residues] of entities.entries()) {
    const isPolymer = residues.size > 1;
    mmcifContent += `${entityId} ${isPolymer ? "polymer" : "non-polymer"} "Chain ${chainId}"\n`;
    entityId++;
  }
  
  // Define atom sites
  mmcifContent += '\n# Atom site information\n';
  mmcifContent += 'loop_\n';
  mmcifContent += '_atom_site.group_PDB\n';
  mmcifContent += '_atom_site.id\n';
  mmcifContent += '_atom_site.type_symbol\n';
  mmcifContent += '_atom_site.label_atom_id\n';
  mmcifContent += '_atom_site.label_alt_id\n';
  mmcifContent += '_atom_site.label_comp_id\n';
  mmcifContent += '_atom_site.label_asym_id\n';
  mmcifContent += '_atom_site.label_entity_id\n';
  mmcifContent += '_atom_site.label_seq_id\n';
  mmcifContent += '_atom_site.pdbx_PDB_ins_code\n';
  mmcifContent += '_atom_site.Cartn_x\n';
  mmcifContent += '_atom_site.Cartn_y\n';
  mmcifContent += '_atom_site.Cartn_z\n';
  mmcifContent += '_atom_site.occupancy\n';
  mmcifContent += '_atom_site.B_iso_or_equiv\n';
  
  // Map chains to entity IDs
  const chainToEntity = new Map();
  let currEntityId = 1;
  for (const chainId of chains) {
    chainToEntity.set(chainId, currEntityId++);
  }
  
  // Add atom site records
  for (const atom of atoms) {
    const entityId = chainToEntity.get(atom.chainId) || 1;
    mmcifContent += `${atom.recordType} ${atom.atomNum} ${atom.element} ${atom.atomName} ${atom.altLoc || '.'} ${atom.residueName} ${atom.chainId} ${entityId} ${atom.residueNum} ${atom.iCode || '.'} ${atom.x.toFixed(3)} ${atom.y.toFixed(3)} ${atom.z.toFixed(3)} ${atom.occupancy.toFixed(2)} ${atom.tempFactor.toFixed(2)}\n`;
  }
  
  return {
    content: mmcifContent,
    atomCount: atoms.length,
    chainCount: chains.size,
    modelCount: models.size || 1,
    note: 'This is a simplified mmCIF conversion focusing on atom coordinates. Some PDB metadata may not be included.'
  };
};

/**
 * Convert mmCIF to PDB format
 * @param {string} content - The mmCIF content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const mmcifToPdb = (content, options) => {
  // Parse mmCIF - this is a very simplified parser
  const lines = content.split('\n');
  
  // Extract structure data blocks
  let title = '';
  const atoms = [];
  
  // Find title
  const titleLine = lines.find(line => line.startsWith('_struct.title'));
  if (titleLine) {
    const titleMatch = titleLine.match(/_struct\.title\s+"(.*)"/);
    if (titleMatch) title = titleMatch[1];
  }
  
  // Find atom_site loop
  let inAtomSiteLoop = false;
  let atomSiteColumnMap = {};
  
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();
    
    if (line === 'loop_') {
      // Check if the next line starts the atom_site loop
      if (i + 1 < lines.length && lines[i + 1].trim().startsWith('_atom_site.')) {
        inAtomSiteLoop = true;
        
        // Parse column definitions
        let columnIndex = 0;
        let j = i + 1;
        
        while (j < lines.length && lines[j].trim().startsWith('_atom_site.')) {
          const columnName = lines[j].trim().substring(11); // Remove '_atom_site.'
          atomSiteColumnMap[columnName] = columnIndex++;
          j++;
        }
        
        // Skip to the data lines
        i = j - 1;
      }
    } else if (inAtomSiteLoop) {
      // Process data lines until we hit another section
      if (line === '' || line.startsWith('_') || line.startsWith('#') || line.startsWith('loop_')) {
        inAtomSiteLoop = false;
      } else {
        // Parse atom data
        const fields = line.split(/\s+/);
        
        if (fields.length >= Object.keys(atomSiteColumnMap).length) {
          const atom = {
            group_PDB: fields[atomSiteColumnMap['group_PDB']] || 'ATOM',
            id: parseInt(fields[atomSiteColumnMap['id']] || '1'),
            type_symbol: fields[atomSiteColumnMap['type_symbol']] || '?',
            label_atom_id: fields[atomSiteColumnMap['label_atom_id']] || '?',
            label_alt_id: fields[atomSiteColumnMap['label_alt_id']] || '.',
            label_comp_id: fields[atomSiteColumnMap['label_comp_id']] || 'UNK',
            label_asym_id: fields[atomSiteColumnMap['label_asym_id']] || 'A',
            label_seq_id: parseInt(fields[atomSiteColumnMap['label_seq_id']] || '1'),
            Cartn_x: parseFloat(fields[atomSiteColumnMap['Cartn_x']] || '0.0'),
            Cartn_y: parseFloat(fields[atomSiteColumnMap['Cartn_y']] || '0.0'),
            Cartn_z: parseFloat(fields[atomSiteColumnMap['Cartn_z']] || '0.0'),
            occupancy: parseFloat(fields[atomSiteColumnMap['occupancy']] || '1.0'),
            B_iso_or_equiv: parseFloat(fields[atomSiteColumnMap['B_iso_or_equiv']] || '0.0')
          };
          
          atoms.push(atom);
        }
      }
    }
  }
  
  // Generate PDB format
  let pdbContent = '';
  
  // Add header
  pdbContent += `HEADER    STRUCTURE FROM MMCIF                     01-JAN-23   NONE\n`;
  
  if (title) {
    pdbContent += `TITLE     ${title}\n`;
  }
  
  // Add atoms
  for (const atom of atoms) {
    const recordType = atom.group_PDB.padEnd(6);
    const atomNum = atom.id.toString().padStart(5);
    const atomName = atom.label_atom_id.padEnd(4);
    const altLoc = atom.label_alt_id === '.' ? ' ' : atom.label_alt_id;
    const residueName = atom.label_comp_id.padEnd(3);
    const chainId = atom.label_asym_id;
    const residueNum = atom.label_seq_id.toString().padStart(4);
    const x = atom.Cartn_x.toFixed(3).padStart(8);
    const y = atom.Cartn_y.toFixed(3).padStart(8);
    const z = atom.Cartn_z.toFixed(3).padStart(8);
    const occupancy = atom.occupancy.toFixed(2).padStart(6);
    const tempFactor = atom.B_iso_or_equiv.toFixed(2).padStart(6);
    const element = atom.type_symbol.padStart(2);
    
    pdbContent += `${recordType}${atomNum} ${atomName}${altLoc}${residueName} ${chainId}${residueNum}    ${x}${y}${z}${occupancy}${tempFactor}          ${element}  \n`;
  }
  
  // Add termination
  pdbContent += 'END\n';
  
  return {
    content: pdbContent,
    atomCount: atoms.length,
    note: 'Simplified conversion from mmCIF to PDB. Some metadata may not be preserved.'
  };
};

/**
 * Convert mmCIF to FASTA (extract protein sequence)
 * @param {string} content - The mmCIF content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const mmcifToFasta = (content, options) => {
  // This is very similar to pdbToFasta, but with different parsing
  const lines = content.split('\n');
  const chains = new Map();
  
  // Amino acid 3-letter to 1-letter conversion
  const aa3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J',
    'XAA': 'X', 'UNK': 'X'
  };
  
  // Extract title
  let title = '';
  const titleLine = lines.find(line => line.startsWith('_struct.title'));
  if (titleLine) {
    const titleMatch = titleLine.match(/_struct\.title\s+"(.*)"/);
    if (titleMatch) title = titleMatch[1];
  }
  
  // Find atom_site loop
  let inAtomSiteLoop = false;
  let atomSiteColumnMap = {};
  
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();
    
    if (line === 'loop_') {
      // Check if the next line starts the atom_site loop
      if (i + 1 < lines.length && lines[i + 1].trim().startsWith('_atom_site.')) {
        inAtomSiteLoop = true;
        
        // Parse column definitions
        let columnIndex = 0;
        let j = i + 1;
        
        while (j < lines.length && lines[j].trim().startsWith('_atom_site.')) {
          const columnName = lines[j].trim().substring(11); // Remove '_atom_site.'
          atomSiteColumnMap[columnName] = columnIndex++;
          j++;
        }
        
        // Skip to the data lines
        i = j - 1;
      }
    } else if (inAtomSiteLoop) {
      // Process data lines until we hit another section
      if (line === '' || line.startsWith('_') || line.startsWith('#') || line.startsWith('loop_')) {
        inAtomSiteLoop = false;
      } else {
        // Parse atom data
        const fields = line.split(/\s+/);
        
        if (fields.length >= Object.keys(atomSiteColumnMap).length) {
          // Check if this is a CA atom (alpha carbon)
          const atomName = fields[atomSiteColumnMap['label_atom_id']] || '';
          const residueName = fields[atomSiteColumnMap['label_comp_id']] || '';
          const chainId = fields[atomSiteColumnMap['label_asym_id']] || '';
          const residueNum = parseInt(fields[atomSiteColumnMap['label_seq_id']] || '0');
          
          if (atomName === 'CA' && residueName in aa3to1) {
            if (!chains.has(chainId)) {
              chains.set(chainId, new Map());
            }
            
            const chain = chains.get(chainId);
            chain.set(residueNum, aa3to1[residueName]);
          }
        }
      }
    }
  }
  
  // Convert to FASTA format
  let fastaContent = '';
  let totalLength = 0;
  
  for (const [chainId, residues] of chains.entries()) {
    // Sort residues by residue number
    const sortedResidues = new Map([...residues.entries()].sort((a, b) => a[0] - b[0]));
    
    // Construct sequence
    const sequence = [...sortedResidues.values()].join('');
    totalLength += sequence.length;
    
    // Create header
    const header = title ? 
      `${title} Chain ${chainId}` : 
      `mmCIF Chain ${chainId}`;
    
    if (options.formatSequence) {
      fastaContent += `>${header}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
    } else {
      fastaContent += `>${header}\n${sequence}\n`;
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount: chains.size,
    totalLength,
    note: 'Sequences were extracted from atom_site records for standard amino acids only.'
  };
};

/**
 * Convert VCF to FASTA (mock implementation)
 * @param {string} content - The VCF content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const vcfToFasta = (content, options) => {
  // This is a placeholder as VCF files don't contain complete sequence data
  // In a real implementation, this would require a reference sequence as additional input
  
  // Parse VCF header
  const lines = content.split('\n');
  let headerLines = [];
  let dataLines = [];
  
  for (const line of lines) {
    if (line.startsWith('#')) {
      headerLines.push(line);
    } else if (line.trim()) {
      dataLines.push(line);
    }
  }
  
  // Extract chromosomes/contigs
  const chroms = new Set();
  for (const line of dataLines) {
    const fields = line.split('\t');
    if (fields.length > 0) {
      chroms.add(fields[0]);
    }
  }
  
  let fastaContent = '';
  
  // Create a placeholder sequence for each chromosome
  for (const chrom of chroms) {
    fastaContent += `>${chrom} (placeholder sequence for VCF variants)\n`;
    
    // Create a dummy sequence of Ns (100 bases)
    const dummySequence = 'N'.repeat(100);
    
    if (options.formatSequence) {
      fastaContent += formatSequenceForFasta(dummySequence, options.lineLength || 60) + '\n';
    } else {
      fastaContent += dummySequence + '\n';
    }
  }
  
  if (chroms.size === 0) {
    fastaContent = '>VCF_derived_sequence (placeholder)\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n';
  }
  
  return {
    content: fastaContent,
    sequenceCount: chroms.size || 1,
    note: 'VCF files contain variant information, not complete sequences. This FASTA contains placeholder sequences only. For a complete conversion, provide the reference sequence file separately.'
  };
};

/**
 * Convert NEWICK to NEXUS format
 * @param {string} content - The NEWICK tree content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const newickToNexus = (content, options) => {
  // Clean up the content
  const newickTree = content.trim();
  
  // Extract taxa from Newick (basic extract - not comprehensive)
  const taxaMatch = newickTree.match(/[A-Za-z0-9_]+/g) || [];
  const uniqueTaxa = [...new Set(taxaMatch)];
  
  // Count number of trees by counting semicolons
  const treeCount = (newickTree.match(/;/g) || []).length;
  
  // Generate NEXUS format
  let nexusContent = '#NEXUS\n\n';
  
  // Add taxa block
  nexusContent += 'BEGIN TAXA;\n';
  nexusContent += `\tDIMENSIONS NTAX=${uniqueTaxa.length};\n`;
  nexusContent += '\tTAXLABELS\n\t\t';
  nexusContent += uniqueTaxa.join('\n\t\t');
  nexusContent += '\n\t;\nEND;\n\n';
  
  // Add trees block
  nexusContent += 'BEGIN TREES;\n';
  
  // Split multiple trees if present
  const trees = newickTree.split(';').filter(tree => tree.trim());
  
  for (let i = 0; i < trees.length; i++) {
    nexusContent += `\tTREE tree${i + 1} = ${trees[i].trim()};\n`;
  }
  
  nexusContent += 'END;\n';
  
  return {
    content: nexusContent,
    taxaCount: uniqueTaxa.length,
    treeCount,
    note: 'Taxa were extracted from the Newick tree, but some complex taxa names with special characters might not be properly identified.'
  };
};

/**
 * Convert NEXUS to NEWICK format
 * @param {string} content - The NEXUS content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const nexusToNewick = (content, options) => {
  // Extract trees from NEXUS file
  const treeBlockMatch = content.match(/BEGIN TREES;([\s\S]*?)END;/i);
  
  if (!treeBlockMatch) {
    throw new Error('No TREES block found in NEXUS file');
  }
  
  const treeBlock = treeBlockMatch[1];
  const treeLines = treeBlock.split('\n')
    .filter(line => line.trim().startsWith('TREE') || line.trim().startsWith('tree'));
  
  if (treeLines.length === 0) {
    throw new Error('No tree definitions found in NEXUS file');
  }
  
  let newickContent = '';
  
  // Extract Newick trees
  for (const line of treeLines) {
    const match = line.match(/(?:TREE|tree)[^=]+=\s*(.*)/i);
    if (match) {
      const tree = match[1].trim();
      // Make sure the tree ends with a semicolon
      newickContent += tree + (tree.endsWith(';') ? '' : ';') + '\n';
    }
  }
  
  // Count number of trees
  const treeCount = (newickContent.match(/;/g) || []).length;
  
  return {
    content: newickContent,
    treeCount,
    note: 'Only tree definitions were extracted. Nexus metadata and other blocks were discarded.'
  };
};

/**
 * Convert SAM to FASTA format
 * @param {string} content - The SAM content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const samToFasta = (content, options) => {
  const lines = content.split('\n');
  let fastaContent = '';
  let sequenceCount = 0;
  let totalLength = 0;
  
  // Skip header lines
  const dataLines = lines.filter(line => !line.startsWith('@') && line.trim());
  
  for (const line of dataLines) {
    const fields = line.split('\t');
    
    if (fields.length >= 10) {
      const readName = fields[0];
      const sequence = fields[9];
      
      if (sequence && sequence !== '*') {
        // Add to FASTA format
        if (options.formatSequence) {
          fastaContent += `>${readName}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
        } else {
          fastaContent += `>${readName}\n${sequence}\n`;
        }
        
        sequenceCount++;
        totalLength += sequence.length;
      }
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount,
    totalLength,
    note: 'Only the sequence field (column 10) was extracted from SAM records. Alignment information was discarded.'
  };
};

/**
 * Convert CLUSTAL to FASTA format
 * @param {string} content - The CLUSTAL content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const clustalToFasta = (content, options) => {
  const lines = content.split('\n');
  const sequences = {};
  let fastaContent = '';
  
  // Find the line with CLUSTAL header
  const headerIndex = lines.findIndex(line => line.includes('CLUSTAL'));
  
  if (headerIndex === -1) {
    throw new Error("Not a valid CLUSTAL format");
  }
  
  // Parse sequences
  for (let i = headerIndex + 1; i < lines.length; i++) {
    const line = lines[i].trim();
    
    // Skip empty lines and conservation lines (starting with spaces)
    if (!line || line.startsWith(' ')) continue;
    
    const parts = line.split(/\s+/);
    
    if (parts.length >= 2) {
      const seqName = parts[0];
      const seqSegment = parts[1];
      
      if (!sequences[seqName]) {
        sequences[seqName] = '';
      }
      
      sequences[seqName] += seqSegment;
    }
  }
  
  // Convert to FASTA
  for (const [name, sequence] of Object.entries(sequences)) {
    if (options.removeGaps) {
      const sequenceWithoutGaps = removeGapsFromSequence(sequence);
      
      if (options.formatSequence) {
        fastaContent += `>${name}\n${formatSequenceForFasta(sequenceWithoutGaps, options.lineLength || 60)}\n`;
      } else {
        fastaContent += `>${name}\n${sequenceWithoutGaps}\n`;
      }
    } else {
      if (options.formatSequence) {
        fastaContent += `>${name}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
      } else {
        fastaContent += `>${name}\n${sequence}\n`;
      }
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount: Object.keys(sequences).length,
    totalLength: Object.values(sequences).reduce((sum, seq) => sum + seq.length, 0),
    note: options.removeGaps ? 'Gaps were removed from the sequences.' : 'Alignment gaps were preserved in the sequences.'
  };
};

/**
 * Convert CLUSTAL to Stockholm format
 * @param {string} content - The CLUSTAL content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const clustalToStockholm = (content, options) => {
  const lines = content.split('\n');
  const sequences = {};
  let conservationLines = [];
  
  // Find the line with CLUSTAL header
  const headerIndex = lines.findIndex(line => line.includes('CLUSTAL'));
  
  if (headerIndex === -1) {
    throw new Error("Not a valid CLUSTAL format");
  }
  
  // Parse sequences and conservation lines
  for (let i = headerIndex + 1; i < lines.length; i++) {
    const line = lines[i].trim();
    
    if (!line) continue;
    
    if (line.startsWith(' ')) {
      // Conservation line
      conservationLines.push(line);
    } else {
      const parts = line.split(/\s+/);
      
      if (parts.length >= 2) {
        const seqName = parts[0];
        const seqSegment = parts[1];
        
        if (!sequences[seqName]) {
          sequences[seqName] = '';
        }
        
        sequences[seqName] += seqSegment;
      }
    }
  }
  
  // Construct Stockholm format
  let stockholmContent = '# STOCKHOLM 1.0\n';
  
  // Add sequences
  for (const [name, sequence] of Object.entries(sequences)) {
    stockholmContent += `${name} ${sequence}\n`;
  }
  
  // Add basic consensus (simplified)
  if (conservationLines.length > 0) {
    // Join conservation symbols
    let consensus = conservationLines.join('').replace(/\s+/g, '');
    
    // Make sure it's the same length as the sequences
    const firstSeqLength = Object.values(sequences)[0]?.length || 0;
    if (consensus.length > firstSeqLength) {
      consensus = consensus.substring(0, firstSeqLength);
    } else if (consensus.length < firstSeqLength) {
      consensus = consensus.padEnd(firstSeqLength, ' ');
    }
    
    stockholmContent += `#=GC cons ${consensus}\n`;
  }
  
  // Add Stockholm end marker
  stockholmContent += '//\n';
  
  return {
    content: stockholmContent,
    sequenceCount: Object.keys(sequences).length,
    alignmentLength: Object.values(sequences)[0]?.length || 0,
    note: 'Converted from CLUSTAL format. Conservation information was preserved where available.'
  };
};

/**
 * Convert Stockholm to FASTA format
 * @param {string} content - The Stockholm content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const stockholmToFasta = (content, options) => {
  const lines = content.split('\n');
  const sequences = {};
  let fastaContent = '';
  
  // Check if it's a valid Stockholm file
  if (!lines.some(line => line.includes('STOCKHOLM'))) {
    throw new Error("Not a valid Stockholm format");
  }
  
  // Parse sequences
  for (const line of lines) {
    // Skip empty lines, comments, and annotations
    if (!line.trim() || line.startsWith('#') || line.startsWith('//')) {
      continue;
    }
    
    // Parse sequence lines
    const parts = line.trim().split(/\s+/);
    
    if (parts.length >= 2) {
      const seqName = parts[0];
      const seqSegment = parts[1];
      
      if (!sequences[seqName]) {
        sequences[seqName] = '';
      }
      
      sequences[seqName] += seqSegment;
    }
  }
  
  // Convert to FASTA
  for (const [name, sequence] of Object.entries(sequences)) {
    if (options.removeGaps) {
      const sequenceWithoutGaps = removeGapsFromSequence(sequence);
      
      if (options.formatSequence) {
        fastaContent += `>${name}\n${formatSequenceForFasta(sequenceWithoutGaps, options.lineLength || 60)}\n`;
      } else {
        fastaContent += `>${name}\n${sequenceWithoutGaps}\n`;
      }
    } else {
      if (options.formatSequence) {
        fastaContent += `>${name}\n${formatSequenceForFasta(sequence, options.lineLength || 60)}\n`;
      } else {
        fastaContent += `>${name}\n${sequence}\n`;
      }
    }
  }
  
  return {
    content: fastaContent,
    sequenceCount: Object.keys(sequences).length,
    totalLength: Object.values(sequences).reduce((sum, seq) => sum + seq.length, 0),
    note: options.removeGaps ? 'Gaps were removed from the sequences.' : 'Alignment gaps were preserved in the sequences.'
  };
};

/**
 * Convert Stockholm to CLUSTAL format
 * @param {string} content - The Stockholm content
 * @param {object} options - Conversion options
 * @returns {object} The conversion result
 */
export const stockholmToClustal = (content, options) => {
  const lines = content.split('\n');
  const sequences = {};
  let consensus = '';
  
  // Check if it's a valid Stockholm file
  if (!lines.some(line => line.includes('STOCKHOLM'))) {
    throw new Error("Not a valid Stockholm format");
  }
  
  // Parse sequences and consensus
  for (const line of lines) {
    if (!line.trim() || line.startsWith('//')) {
      continue;
    }
    
    if (line.startsWith('#=GC')) {
      // Consensus line
      const parts = line.split(/\s+/);
      if (parts.length >= 3 && parts[1] === 'cons') {
        consensus = parts[2];
      }
    } else if (!line.startsWith('#')) {
      // Sequence line
      const parts = line.trim().split(/\s+/);
      
      if (parts.length >= 2) {
        const seqName = parts[0];
        const seqSegment = parts[1];
        
        if (!sequences[seqName]) {
          sequences[seqName] = '';
        }
