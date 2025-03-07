// FormatDetection.js - Handles detection and validation of biological file formats

/**
 * Detect the format of a biological data file based on its content and name
 * @param {string} content - The content of the file
 * @param {string} fileName - The name of the file
 * @param {object} formatDefinitions - The available format definitions
 * @returns {string} The detected format key, or null if not detected
 */
export const detectFileFormat = (content, fileName, formatDefinitions) => {
  // Early return for empty content
  if (!content || content.trim() === '') {
    return null;
  }
  
  // First, check file extension
  const extension = fileName ? ('.' + fileName.split('.').pop().toLowerCase()) : '';
  
  // Find formats matching the extension
  const formatsByExtension = [];
  
  for (const [key, format] of Object.entries(formatDefinitions)) {
    if (format.extensions && format.extensions.includes(extension)) {
      formatsByExtension.push(key);
    }
  }
  
  // If only one format matches by extension, verify content
  if (formatsByExtension.length === 1) {
    const format = formatsByExtension[0];
    if (matchesFormatContent(content, format)) {
      return format;
    }
  }
  
  // Try to detect by content
  return detectFormatByContent(content);
};

/**
 * Detect the format based solely on content
 * @param {string} content - The content of the file
 * @returns {string} The detected format key, or null if not detected
 */
export const detectFormatByContent = (content) => {
  const trimmedContent = content.trim();
  
  // FASTA format (starts with >)
  if (trimmedContent.startsWith('>')) {
    return 'fasta';
  }
  
  // FASTQ format (starts with @ and has 4-line pattern)
  if (trimmedContent.startsWith('@')) {
    const lines = trimmedContent.split('\n');
    if (lines.length >= 4 && lines[2].startsWith('+')) {
      return 'fastq';
    }
  }
  
  // GenBank format (has LOCUS at start)
  if (trimmedContent.startsWith('LOCUS ')) {
    return 'genbank';
  }
  
  // EMBL format (starts with ID or has EMBL format markers)
  if (trimmedContent.startsWith('ID ') && trimmedContent.includes('SQ ')) {
    return 'embl';
  }
  
  // PDB format (has ATOM or HETATM records)
  if (trimmedContent.includes('ATOM  ') || trimmedContent.includes('HETATM')) {
    return 'pdb';
  }
  
  // mmCIF format (has data_ and loop_)
  if (trimmedContent.startsWith('data_') && trimmedContent.includes('loop_')) {
    return 'mmcif';
  }
  
  // GFF format (has ##gff-version)
  if (trimmedContent.includes('##gff-version')) {
    return 'gff';
  }
  
  // VCF format (has ##fileformat=VCF)
  if (trimmedContent.includes('##fileformat=VCF')) {
    return 'vcf';
  }
  
  // BED format (tab-delimited with chromosome, start, end)
  if (/^[^\s]+\t\d+\t\d+/.test(trimmedContent)) {
    return 'bed';
  }
  
  // SAM format (has @HD header)
  if (trimmedContent.includes('@HD\tVN:')) {
    return 'sam';
  }
  
  // CLUSTAL format
  if (trimmedContent.includes('CLUSTAL W') || trimmedContent.includes('CLUSTAL Multiple Alignment')) {
    return 'clustal';
  }
  
  // Stockholm format
  if (trimmedContent.includes('# STOCKHOLM') && trimmedContent.includes('//')) {
    return 'stockholm';
  }
  
  // Newick format (tree format with parentheses and semicolon)
  if (/^\s*\(.*\);\s*$/.test(trimmedContent)) {
    return 'newick';
  }
  
  // NEXUS format (has #NEXUS)
  if (trimmedContent.startsWith('#NEXUS')) {
    return 'nexus';
  }
  
  // Could not determine format
  return null;
};

/**
 * Check if content matches a specific format's expected patterns
 * @param {string} content - The content to check
 * @param {string} format - The format key to check against
 * @returns {boolean} Whether the content matches the format
 */
export const matchesFormatContent = (content, format) => {
  if (!content) return false;
  
  const trimmedContent = content.trim();
  
  switch (format) {
    case 'fasta':
      return trimmedContent.startsWith('>');
      
    case 'fastq':
      const fastqLines = trimmedContent.split('\n');
      return fastqLines.length >= 4 && 
             fastqLines[0].startsWith('@') && 
             fastqLines[2].startsWith('+');
      
    case 'genbank':
      return trimmedContent.startsWith('LOCUS ') || 
             trimmedContent.includes('FEATURES') && 
             trimmedContent.includes('ORIGIN');
      
    case 'embl':
      return trimmedContent.startsWith('ID ') && 
             trimmedContent.includes('SQ ');
      
    case 'pdb':
      return trimmedContent.includes('ATOM  ') || 
             trimmedContent.includes('HETATM');
      
    case 'mmcif':
      return trimmedContent.startsWith('data_') && 
             trimmedContent.includes('loop_');
      
    case 'gff':
      return trimmedContent.includes('##gff-version') || 
             /^[^\s#][^\t]*\t[^\t]*\t[^\t]*\t\d+\t\d+/.test(trimmedContent);
      
    case 'bed':
      return /^[^\s]+\t\d+\t\d+/.test(trimmedContent);
      
    case 'vcf':
      return trimmedContent.includes('##fileformat=VCF');
      
    case 'sam':
      return trimmedContent.includes('@HD\tVN:') || 
             trimmedContent.includes('@SQ\tSN:');
      
    case 'clustal':
      return trimmedContent.includes('CLUSTAL W') || 
             trimmedContent.includes('CLUSTAL Multiple Alignment');
      
    case 'stockholm':
      return trimmedContent.includes('# STOCKHOLM') && 
             trimmedContent.includes('//');
      
    case 'newick':
      return /^\s*\(.*\);\s*$/.test(trimmedContent);
      
    case 'nexus':
      return trimmedContent.startsWith('#NEXUS');
      
    default:
      return false;
  }
};

/**
 * Validate if a file's content conforms to a specific format
 * @param {string} content - The content to validate
 * @param {string} format - The format to validate against
 * @returns {object} Validation result with valid flag and any error messages
 */
export const validateFormat = (content, format) => {
  if (!content) {
    return { valid: false, errors: ['Empty content'] };
  }
  
  const errors = [];
  
  switch (format) {
    case 'fasta':
      if (!content.startsWith('>')) {
        errors.push('FASTA file must start with ">" character');
      }
      
      // Check if all sequences have headers
      const fastaLines = content.split('\n');
      let hasHeader = false;
      
      for (let i = 0; i < fastaLines.length; i++) {
        const line = fastaLines[i].trim();
        
        if (line.startsWith('>')) {
          hasHeader = true;
        } else if (line && !hasHeader) {
          errors.push('Found sequence data without a header');
          break;
        }
      }
      break;
      
    case 'fastq':
      const fastqLines = content.split('\n');
      
      // Check if the file follows the 4-line pattern (header, sequence, +, quality)
      for (let i = 0; i < fastqLines.length; i += 4) {
        if (i + 3 >= fastqLines.length) {
          if (i < fastqLines.length && fastqLines[i].trim()) {
            errors.push('Incomplete FASTQ record at the end of the file');
          }
          break;
        }
        
        if (!fastqLines[i].startsWith('@')) {
          errors.push(`Line ${i+1} should start with '@'`);
        }
        
        if (!fastqLines[i+2].startsWith('+')) {
          errors.push(`Line ${i+3} should start with '+'`);
        }
        
        // Check if sequence length matches quality length
        const seqLength = fastqLines[i+1].trim().length;
        const qualLength = fastqLines[i+3].trim().length;
        
        if (seqLength !== qualLength) {
          errors.push(`Sequence and quality length mismatch at line ${i+1}`);
        }
      }
      break;
      
    case 'genbank':
      if (!content.includes('LOCUS ')) {
        errors.push('GenBank file must contain a LOCUS line');
      }
      
      if (!content.includes('ORIGIN')) {
        errors.push('GenBank file must contain an ORIGIN section');
      }
      
      if (!content.includes('//')) {
        errors.push('GenBank file must end with "//"');
      }
      break;
      
    case 'embl':
      if (!content.startsWith('ID ')) {
        errors.push('EMBL file must start with an ID line');
      }
      
      if (!content.includes('SQ ')) {
        errors.push('EMBL file must contain an SQ section');
      }
      
      if (!content.includes('//')) {
        errors.push('EMBL file must end with "//"');
      }
      break;
      
    case 'pdb':
      if (!content.includes('ATOM  ') && !content.includes('HETATM')) {
        errors.push('PDB file must contain ATOM or HETATM records');
      }
      
      // Check for END marker
      if (!content.includes('END')) {
        errors.push('PDB file should end with an END record');
      }
      break;
      
    case 'mmcif':
      if (!content.startsWith('data_')) {
        errors.push('mmCIF file must start with a data_ block');
      }
      
      if (!content.includes('loop_')) {
        errors.push('mmCIF file must contain at least one loop_ section');
      }
      break;
      
    case 'gff':
      if (!content.includes('##gff-version')) {
        errors.push('GFF file should include a ##gff-version pragma');
      }
      
      // Check a sample of data lines
      const gffLines = content.split('\n').filter(line => !line.startsWith('#') && line.trim());
      
      for (let i = 0; i < Math.min(gffLines.length, 5); i++) {
        const parts = gffLines[i].split('\t');
        if (parts.length < 8) {
          errors.push(`GFF line does not have the required 8+ columns: ${gffLines[i]}`);
          break;
        }
      }
      break;
      
    case 'bed':
      // Check a sample of data lines
      const bedLines = content.split('\n').filter(line => !line.startsWith('#') && line.trim());
      
      for (let i = 0; i < Math.min(bedLines.length, 5); i++) {
        const parts = bedLines[i].split('\t');
        if (parts.length < 3) {
          errors.push(`BED line does not have the required 3+ columns: ${bedLines[i]}`);
          break;
        }
        
        if (isNaN(parseInt(parts[1])) || isNaN(parseInt(parts[2]))) {
          errors.push(`BED line has invalid start/end positions: ${bedLines[i]}`);
          break;
        }
      }
      break;
      
    case 'vcf':
      if (!content.includes('##fileformat=VCF')) {
        errors.push('VCF file must include a ##fileformat=VCF header');
      }
      
      if (!content.includes('#CHROM\tPOS\tID\tREF\tALT')) {
        errors.push('VCF file must include the standard column header line');
      }
      break;
      
    case 'sam':
      // Check for SAM header
      if (!content.includes('@HD') && !content.includes('@SQ')) {
        errors.push('SAM file should include @HD or @SQ header lines');
      }
      
      // Check a sample of data lines
      const samLines = content.split('\n').filter(line => !line.startsWith('@') && line.trim());
      
      for (let i = 0; i < Math.min(samLines.length, 5); i++) {
        const parts = samLines[i].split('\t');
        if (parts.length < 11) {
          errors.push(`SAM alignment line does not have the required 11+ columns: ${samLines[i]}`);
          break;
        }
      }
      break;
      
    case 'clustal':
      if (!content.includes('CLUSTAL W') && !content.includes('CLUSTAL Multiple Alignment')) {
        errors.push('CLUSTAL file must include a CLUSTAL header line');
      }
      break;
      
    case 'stockholm':
      if (!content.includes('# STOCKHOLM')) {
        errors.push('Stockholm file must include a # STOCKHOLM header');
      }
      
      if (!content.includes('//')) {
        errors.push('Stockholm file must end with "//"');
      }
      break;
      
    case 'newick':
      if (!content.includes('(') || !content.includes(')') || !content.includes(';')) {
        errors.push('Newick tree format must include parentheses and end with a semicolon');
      }
      
      // Basic check for balanced parentheses
      let count = 0;
      for (const char of content) {
        if (char === '(') count++;
        if (char === ')') count--;
        if (count < 0) {
          errors.push('Unbalanced parentheses in Newick tree');
          break;
        }
      }
      
      if (count !== 0) {
        errors.push('Unbalanced parentheses in Newick tree');
      }
      break;
      
    case 'nexus':
      if (!content.startsWith('#NEXUS')) {
        errors.push('NEXUS file must start with #NEXUS');
      }
      
      if (!content.includes('BEGIN') || !content.includes('END;')) {
        errors.push('NEXUS file must include at least one BEGIN/END block');
      }
      break;
      
    default:
      errors.push(`Unknown format: ${format}`);
      break;
  }
  
  return {
    valid: errors.length === 0,
    errors: errors
  };
};

/**
 * Detect if a sequence is likely DNA, RNA, or protein
 * @param {string} sequence - The sequence to analyze
 * @returns {string} The detected sequence type ('dna', 'rna', 'protein', or 'unknown')
 */
export const detectSequenceType = (sequence) => {
  if (!sequence || typeof sequence !== 'string') {
    return 'unknown';
  }
  
  const upperSeq = sequence.toUpperCase();
  
  // Count bases/residues
  const counts = {};
  for (const char of upperSeq) {
    if (!/[A-Z]/.test(char)) continue; // Skip non-alphabetic characters
    counts[char] = (counts[char] || 0) + 1;
  }
  
  const totalChars = Object.values(counts).reduce((sum, count) => sum + count, 0);
  if (totalChars === 0) return 'unknown';
  
  // Calculate percentage of DNA/RNA bases
  const dnaRnaBases = ['A', 'C', 'G', 'T', 'U', 'N'];
  const dnaRnaBaseCount = dnaRnaBases.reduce((sum, base) => sum + (counts[base] || 0), 0);
  const dnaRnaPercent = (dnaRnaBaseCount / totalChars) * 100;
  
  // Check for T and U to differentiate DNA and RNA
  const tPercent = ((counts['T'] || 0) / totalChars) * 100;
  const uPercent = ((counts['U'] || 0) / totalChars) * 100;
  
  if (dnaRnaPercent > 90) {
    // Likely nucleotide sequence
    if (tPercent > 1 && uPercent < 1) {
      return 'dna';
    } else if (uPercent > 1 && tPercent < 1) {
      return 'rna';
    } else {
      // Default to DNA if unclear
      return 'dna';
    }
  } else {
    // Check for common amino acids that aren't in DNA/RNA
    const proteinSpecificAA = ['E', 'F', 'I', 'L', 'P', 'Q', 'Y'];
    const proteinSpecificCount = proteinSpecificAA.reduce((sum, aa) => sum + (counts[aa] || 0), 0);
    
    if (proteinSpecificCount > 0) {
      return 'protein';
    }
  }
  
  return 'unknown';
};
