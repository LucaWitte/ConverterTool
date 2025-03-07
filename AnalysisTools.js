// AnalysisTools.js - Provides analysis tools for different biological data formats

import { extractFastaSequences, normalizeNucleotideSequence, translateToProtein } from './FormatConverters';
import { detectSequenceType } from './FormatDetection';

/**
 * Available analysis tools by format
 */
export const AnalysisTools = {
  // FASTA format tools
  fasta: [
    {
      id: 'sequence_stats',
      name: 'Sequence Statistics',
      description: 'Basic statistics about sequences'
    },
    {
      id: 'gc_content',
      name: 'GC Content Analysis',
      description: 'Calculate GC content for each sequence'
    },
    {
      id: 'translate',
      name: 'Translate to Protein',
      description: 'Translate nucleotide sequences to protein'
    },
    {
      id: 'find_motifs',
      name: 'Find Sequence Motifs',
      description: 'Search for specific sequence patterns'
    }
  ],
  
  // FASTQ format tools
  fastq: [
    {
      id: 'quality_stats',
      name: 'Quality Statistics',
      description: 'Analyze quality scores distribution'
    },
    {
      id: 'sequence_stats',
      name: 'Sequence Statistics',
      description: 'Basic statistics about sequences'
    },
    {
      id: 'gc_content',
      name: 'GC Content Analysis',
      description: 'Calculate GC content for each sequence'
    }
  ],
  
  // GenBank format tools
  genbank: [
    {
      id: 'feature_stats',
      name: 'Feature Statistics',
      description: 'Analyze feature types and distribution'
    },
    {
      id: 'sequence_stats',
      name: 'Sequence Statistics',
      description: 'Basic statistics about sequences'
    },
    {
      id: 'gc_content',
      name: 'GC Content Analysis',
      description: 'Calculate GC content for regions'
    }
  ],
  
  // PDB format tools
  pdb: [
    {
      id: 'structure_analysis',
      name: 'Structure Analysis',
      description: 'Analyze protein structure properties'
    },
    {
      id: 'residue_distribution',
      name: 'Residue Distribution',
      description: 'Analyze amino acid frequencies'
    }
  ],
  
  // CLUSTAL format tools
  clustal: [
    {
      id: 'conservation',
      name: 'Conservation Analysis',
      description: 'Analyze sequence conservation'
    },
    {
      id: 'sequence_stats',
      name: 'Alignment Statistics',
      description: 'Basic statistics about the alignment'
    }
  ],
  
  // Other formats
  embl: [
    {
      id: 'feature_stats',
      name: 'Feature Statistics',
      description: 'Analyze feature types and distribution'
    },
    {
      id: 'sequence_stats',
      name: 'Sequence Statistics',
      description: 'Basic statistics about sequences'
    }
  ],
  
  stockholm: [
    {
      id: 'conservation',
      name: 'Conservation Analysis',
      description: 'Analyze sequence conservation'
    }
  ]
};

/**
 * Main function to perform analysis based on the selected tool
 * @param {string} toolId - ID of the selected analysis tool
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @param {function} progressCallback - Callback for progress updates
 * @param {function} successCallback - Callback for successful analysis
 * @param {function} errorCallback - Callback for error handling
 */
export const performAnalysis = (
  toolId,
  content,
  format,
  progressCallback,
  successCallback,
  errorCallback
) => {
  try {
    // Update progress
    progressCallback(10);
    
    // Map of analysis functions
    const analysisFunctions = {
      'sequence_stats': analyzeSequenceStats,
      'gc_content': analyzeGCContent,
      'translate': translateNucleotide,
      'find_motifs': findMotifs,
      'quality_stats': analyzeQualityStats,
      'feature_stats': analyzeFeatureStats,
      'structure_analysis': analyzeStructure,
      'residue_distribution': analyzeResidueDistribution,
      'conservation': analyzeConservation
    };
    
    // Check if analysis is supported
    if (!analysisFunctions[toolId]) {
      throw new Error(`Analysis tool ${toolId} is not supported`);
    }
    
    // Update progress
    progressCallback(30);
    
    // Perform the analysis
    const result = analysisFunctions[toolId](content, format);
    
    // Update progress
    progressCallback(80);
    
    // Return the result
    successCallback({
      type: toolId,
      title: getToolTitle(toolId, format),
      data: result
    });
    
  } catch (error) {
    errorCallback(`Analysis error: ${error.message}`);
  }
};

/**
 * Get the title of an analysis tool
 * @param {string} toolId - ID of the analysis tool
 * @param {string} format - Format of the content
 * @returns {string} Tool title
 */
const getToolTitle = (toolId, format) => {
  const tools = AnalysisTools[format] || [];
  const tool = tools.find(t => t.id === toolId);
  return tool ? tool.name : toolId;
};

/**
 * Analyze basic sequence statistics
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeSequenceStats = (content, format) => {
  let sequences = [];
  let sequenceType = 'unknown';
  
  if (format === 'fasta') {
    sequences = extractFastaSequences(content);
    if (sequences.length > 0) {
      sequenceType = detectSequenceType(sequences[0].sequence);
    }
  } else if (format === 'fastq') {
    // Extract sequences from FASTQ
    const lines = content.split('\n');
    for (let i = 0; i < lines.length; i += 4) {
      if (i + 3 < lines.length && lines[i].startsWith('@')) {
        const header = lines[i].substring(1).trim();
        const sequence = lines[i + 1].trim();
        
        if (sequence) {
          sequences.push({
            header,
            sequence,
            length: sequence.length
          });
        }
      }
    }
    
    if (sequences.length > 0) {
      sequenceType = detectSequenceType(sequences[0].sequence);
    }
  } else if (format === 'genbank') {
    // Extract sequence from GenBank
    const sections = content.split('//');
    
    for (const section of sections) {
      if (!section.trim()) continue;
      
      // Extract LOCUS information
      const locusMatch = section.match(/LOCUS\s+(\S+)/);
      const locus = locusMatch ? locusMatch[1] : 'Unknown';
      
      // Extract sequence
      const originMatch = section.match(/ORIGIN([\s\S]*?)(?=$)/);
      if (originMatch) {
        let sequence = originMatch[1].replace(/\d+|\s+/g, '');
        
        if (sequence) {
          sequences.push({
            header: locus,
            sequence,
            length: sequence.length
          });
        }
      }
    }
    
    if (sequences.length > 0) {
      sequenceType = detectSequenceType(sequences[0].sequence);
    }
  } else if (format === 'clustal') {
    // Extract sequences from CLUSTAL
    const lines = content.split('\n');
    const seqMap = {};
    
    // Find the line with CLUSTAL header
    const headerIndex = lines.findIndex(line => line.includes('CLUSTAL'));
    
    if (headerIndex !== -1) {
      for (let i = headerIndex + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        
        // Skip empty lines and conservation lines (starting with spaces)
        if (!line || line.startsWith(' ')) continue;
        
        const parts = line.split(/\s+/);
        
        if (parts.length >= 2) {
          const seqName = parts[0];
          const seqSegment = parts[1];
          
          if (!seqMap[seqName]) {
            seqMap[seqName] = '';
          }
          
          seqMap[seqName] += seqSegment;
        }
      }
      
      // Convert to array
      for (const [name, sequence] of Object.entries(seqMap)) {
        sequences.push({
          header: name,
          sequence,
          length: sequence.length
        });
      }
      
      if (sequences.length > 0) {
        sequenceType = detectSequenceType(sequences[0].sequence);
      }
    }
  }
  
  // Calculate statistics
  const lengths = sequences.map(seq => seq.length);
  const totalLength = lengths.reduce((sum, len) => sum + len, 0);
  const minLength = Math.min(...lengths);
  const maxLength = Math.max(...lengths);
  const averageLength = sequences.length > 0 ? totalLength / sequences.length : 0;
  
  // Calculate base/residue composition
  const composition = {};
  let totalBases = 0;
  
  for (const seq of sequences) {
    for (const base of seq.sequence.toUpperCase()) {
      if (/[A-Z]/.test(base)) {
        composition[base] = (composition[base] || 0) + 1;
        totalBases++;
      }
    }
  }
  
  // Calculate percentage composition
  const compositionPercent = {};
  for (const [base, count] of Object.entries(composition)) {
    compositionPercent[base] = `${((count / totalBases) * 100).toFixed(2)}%`;
  }
  
  // Calculate GC content for nucleotide sequences
  let gcContent = null;
  if (sequenceType === 'dna' || sequenceType === 'rna') {
    const gc = (composition['G'] || 0) + (composition['C'] || 0);
    const gcPercent = (gc / totalBases) * 100;
    gcContent = `${gcPercent.toFixed(2)}%`;
  }
  
  // Prepare sequence previews
  const seqPreviews = sequences.slice(0, 50).map(seq => ({
    header: seq.header,
    length: seq.length,
    preview: seq.sequence.length > 50 ? 
      seq.sequence.substring(0, 25) + '...' + seq.sequence.substring(seq.sequence.length - 25) : 
      seq.sequence
  }));
  
  return {
    sequenceType,
    stats: {
      sequenceCount: sequences.length,
      totalLength,
      minLength,
      maxLength,
      averageLength: averageLength.toFixed(1)
    },
    composition,
    compositionPercent,
    gcContent,
    sequences: seqPreviews
  };
};

/**
 * Analyze GC content
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeGCContent = (content, format) => {
  let sequences = [];
  
  if (format === 'fasta') {
    sequences = extractFastaSequences(content);
  } else if (format === 'fastq') {
    // Extract sequences from FASTQ
    const lines = content.split('\n');
    for (let i = 0; i < lines.length; i += 4) {
      if (i + 3 < lines.length && lines[i].startsWith('@')) {
        const header = lines[i].substring(1).trim();
        const sequence = lines[i + 1].trim();
        
        if (sequence) {
          sequences.push({
            header,
            sequence,
            length: sequence.length
          });
        }
      }
    }
  } else if (format === 'genbank') {
    // Extract sequence from GenBank
    const sections = content.split('//');
    
    for (const section of sections) {
      if (!section.trim()) continue;
      
      // Extract LOCUS information
      const locusMatch = section.match(/LOCUS\s+(\S+)/);
      const locus = locusMatch ? locusMatch[1] : 'Unknown';
      
      // Extract sequence
      const originMatch = section.match(/ORIGIN([\s\S]*?)(?=$)/);
      if (originMatch) {
        let sequence = originMatch[1].replace(/\d+|\s+/g, '');
        
        if (sequence) {
          sequences.push({
            header: locus,
            sequence,
            length: sequence.length
          });
        }
      }
    }
  }
  
  // Calculate GC content for each sequence
  const sequenceGcContents = [];
  let totalGc = 0;
  let totalBases = 0;
  
  for (const seq of sequences) {
    const upperSeq = seq.sequence.toUpperCase();
    let gcCount = 0;
    let validBases = 0;
    
    for (const base of upperSeq) {
      if (base === 'G' || base === 'C') {
        gcCount++;
        validBases++;
      } else if (base === 'A' || base === 'T' || base === 'U') {
        validBases++;
      }
    }
    
    const gcPercent = validBases > 0 ? (gcCount / validBases) * 100 : 0;
    
    sequenceGcContents.push({
      header: seq.header,
      length: seq.length,
      gcCount,
      gcContent: `${gcPercent.toFixed(2)}%`
    });
    
    totalGc += gcCount;
    totalBases += validBases;
  }
  
  // Calculate overall GC content
  const overallGcPercent = totalBases > 0 ? (totalGc / totalBases) * 100 : 0;
  const overallGcContent = `${overallGcPercent.toFixed(2)}%`;
  
  return {
    sequenceCount: sequences.length,
    overallGcContent,
    sequenceGcContents: sequenceGcContents.sort((a, b) => {
      const gcA = parseFloat(a.gcContent);
      const gcB = parseFloat(b.gcContent);
      return gcB - gcA; // Sort by GC content descending
    }),
    note: 'GC content calculated as the percentage of G and C bases relative to all valid bases (A, T, G, C, U)'
  };
};

/**
 * Translate nucleotide sequences to protein
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const translateNucleotide = (content, format) => {
  if (format !== 'fasta') {
    throw new Error('Translation is only supported for FASTA format');
  }
  
  const sequences = extractFastaSequences(content);
  const translatedSequences = [];
  let fastaFormat = '';
  
  for (const seq of sequences) {
    // Normalize the sequence (uppercase and DNA)
    const normalizedSeq = normalizeNucleotideSequence(seq.sequence);
    
    // Translate to protein
    const protein = translateToProtein(normalizedSeq);
    
    translatedSequences.push({
      header: seq.header,
      nucleotideLength: normalizedSeq.length,
      proteinLength: protein.length
    });
    
    // Add to FASTA format
    fastaFormat += `>${seq.header} [Translated]\n`;
    
    // Format protein sequence with 60 chars per line
    for (let i = 0; i < protein.length; i += 60) {
      fastaFormat += protein.substring(i, i + 60) + '\n';
    }
  }
  
  return {
    sequenceCount: sequences.length,
    translatedSequences,
    fastaFormat
  };
};

/**
 * Search for sequence motifs
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @param {string} motif - Motif to search for (option)
 * @returns {object} Analysis results
 */
const findMotifs = (content, format, motif = null) => {
  if (format !== 'fasta') {
    throw new Error('Motif finding is only supported for FASTA format');
  }
  
  const sequences = extractFastaSequences(content);
  
  // Common motifs to search for if none specified
  const motifs = motif ? [motif] : [
    // DNA motifs
    'TATA',         // TATA box
    'CAAT',         // CAAT box
    'AATAAA',       // Polyadenylation signal
    'GGTAAG',       // Splice donor
    'TTTTTTTT',     // Poly-T
    // Protein motifs
    'KR',           // Basic residues
    'DE',           // Acidic residues
    'DKTGT',        // P-type ATPase phosphorylation site
    'N[^P][ST][^P]' // N-glycosylation site
  ];
  
  const motifResults = [];
  
  for (const currentMotif of motifs) {
    const results = {
      motif: currentMotif,
      totalOccurrences: 0,
      sequenceHits: []
    };
    
    // Create a regex for the motif
    let motifRegex;
    try {
      motifRegex = new RegExp(currentMotif, 'gi');
    } catch (e) {
      // If regex creation fails, use literal string
      motifRegex = new RegExp(currentMotif.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&'), 'gi');
    }
    
    for (const seq of sequences) {
      // Find all matches
      const matches = [...seq.sequence.matchAll(motifRegex)];
      
      if (matches.length > 0) {
        results.totalOccurrences += matches.length;
        
        results.sequenceHits.push({
          header: seq.header,
          occurrences: matches.length,
          positions: matches.map(match => match.index + 1) // 1-based positions
        });
      }
    }
    
    if (results.totalOccurrences > 0) {
      motifResults.push(results);
    }
  }
  
  return {
    sequenceCount: sequences.length,
    motifResults: motifResults.sort((a, b) => b.totalOccurrences - a.totalOccurrences)
  };
};

/**
 * Analyze quality statistics from FASTQ
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeQualityStats = (content, format) => {
  if (format !== 'fastq') {
    throw new Error('Quality analysis is only supported for FASTQ format');
  }
  
  const lines = content.split('\n');
  const reads = [];
  
  // Extract quality scores
  for (let i = 0; i < lines.length; i += 4) {
    if (i + 3 < lines.length && lines[i].startsWith('@')) {
      const header = lines[i].substring(1).trim();
      const sequence = lines[i + 1].trim();
      const quality = lines[i + 3].trim();
      
      if (sequence && quality && sequence.length === quality.length) {
        reads.push({
          header,
          sequence,
          quality,
          length: sequence.length
        });
      }
    }
  }
  
  // Calculate quality statistics
  const positionQuality = {};
  const qualityDistribution = {};
  let totalQualityScore = 0;
  let totalBases = 0;
  let minQuality = 999;
  let maxQuality = 0;
  
  // Process each read
  for (const read of reads) {
    for (let i = 0; i < read.quality.length; i++) {
      const qualChar = read.quality[i];
      const qualScore = qualChar.charCodeAt(0) - 33; // Convert from ASCII to Phred
      
      // Update position-specific quality
      if (!positionQuality[i]) {
        positionQuality[i] = {
          sum: 0,
          count: 0,
          min: 999,
          max: 0
        };
      }
      
      positionQuality[i].sum += qualScore;
      positionQuality[i].count++;
      positionQuality[i].min = Math.min(positionQuality[i].min, qualScore);
      positionQuality[i].max = Math.max(positionQuality[i].max, qualScore);
      
      // Update quality distribution
      qualityDistribution[qualScore] = (qualityDistribution[qualScore] || 0) + 1;
      
      // Update overall stats
      totalQualityScore += qualScore;
      totalBases++;
      minQuality = Math.min(minQuality, qualScore);
      maxQuality = Math.max(maxQuality, qualScore);
    }
  }
  
  // Calculate average quality by position
  const positionQualityArray = Object.entries(positionQuality).map(([pos, stats]) => ({
    position: parseInt(pos) + 1, // 1-based position
    mean: stats.sum / stats.count,
    min: stats.min,
    max: stats.max
  }));
  
  // Calculate quality distribution percentages
  const qualityDistributionArray = Object.entries(qualityDistribution).map(([score, count]) => ({
    score: parseInt(score),
    count,
    percentage: (count / totalBases) * 100
  })).sort((a, b) => a.score - b.score);
  
  // Calculate read length distribution
  const readLengths = {};
  for (const read of reads) {
    readLengths[read.length] = (readLengths[read.length] || 0) + 1;
  }
  
  const readLengthDistribution = Object.entries(readLengths).map(([length, count]) => ({
    length: parseInt(length),
    count,
    percentage: (count / reads.length) * 100
  })).sort((a, b) => a.length - b.length);
  
  return {
    readCount: reads.length,
    overallStats: {
      meanQuality: totalBases > 0 ? (totalQualityScore / totalBases).toFixed(2) : 0,
      minQuality,
      maxQuality,
      totalBases
    },
    positionQuality: positionQualityArray,
    qualityDistribution: qualityDistributionArray,
    readLengthDistribution
  };
};

/**
 * Analyze feature statistics from GenBank or EMBL
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeFeatureStats = (content, format) => {
  if (format !== 'genbank' && format !== 'embl') {
    throw new Error('Feature analysis is only supported for GenBank and EMBL formats');
  }
  
  const features = [];
  let metadata = {};
  
  if (format === 'genbank') {
    // Extract basic metadata
    const locusMatch = content.match(/LOCUS\s+(\S+)\s+(\d+)\s+bp/);
    if (locusMatch) {
      metadata.name = locusMatch[1];
      metadata.length = parseInt(locusMatch[2]);
    }
    
    const definitionMatch = content.match(/DEFINITION\s+(.*?)(?=\n\S)/s);
    if (definitionMatch) {
      metadata.definition = definitionMatch[1].replace(/\n\s+/g, ' ').trim();
    }
    
    const organismMatch = content.match(/\s+ORGANISM\s+(.*?)(?=\n\S)/s);
    if (organismMatch) {
      metadata.organism = organismMatch[1].replace(/\n\s+/g, ' ').trim();
    }
    
    // Extract features
    const featuresMatch = content.match(/FEATURES\s+Location\/Qualifiers\s+(.*?)(?=ORIGIN|\n\S\S)/s);
    if (featuresMatch) {
      const featureText = featuresMatch[1];
      const featureLines = featureText.split('\n');
      
      let currentFeature = null;
      
      for (const line of featureLines) {
        if (line.match(/^\s{5}\S/)) {
          // New feature
          if (currentFeature) {
            features.push(currentFeature);
          }
          
          const parts = line.trim().split(/\s+/, 2);
          const type = parts[0];
          const location = parts.length > 1 ? parts[1] : '';
          currentFeature = { type, location, qualifiers: {} };
        } else if (line.match(/^\s{21}\//) && currentFeature) {
          // Qualifier
          const match = line.trim().match(/\/(\S+)=(?:"([^"]*)"|([\S]+))/);
          if (match) {
            const key = match[1];
            const value = match[2] || match[3];
            currentFeature.qualifiers[key] = value;
          }
        }
      }
      
      if (currentFeature) {
        features.push(currentFeature);
      }
    }
  } else if (format === 'embl') {
    // Extract ID information
    const idMatch = content.match(/ID\s+([^;]+);\s+SV\s+(\d+);\s+.*?;\s+(\d+)\s+BP/);
    if (idMatch) {
      metadata.name = idMatch[1].trim();
      metadata.version = idMatch[2];
      metadata.length = parseInt(idMatch[3]);
    }
    
    // Extract DE (description)
    const deMatch = content.match(/DE\s+(.*?)(?=\nXX|\n\S\S)/s);
    if (deMatch) {
      metadata.definition = deMatch[1].replace(/\n\s+/g, ' ').trim();
    }
    
    // Extract OS (organism)
    const osMatch = content.match(/OS\s+(.*?)(?=\nXX|\n\S\S)/s);
    if (osMatch) {
      metadata.organism = osMatch[1].replace(/\n\s+/g, ' ').trim();
    }
    
    // Extract features
    const featureBlock = content.match(/FH.*?\nFT([\s\S]*?)(?=XX|$)/s);
    
    if (featureBlock) {
      const featureLines = featureBlock[1].split('\n');
      
      let currentFeature = null;
      
      for (let i = 0; i < featureLines.length; i++) {
        const line = featureLines[i];
        
        if (line.match(/^\s+\S+\s+\S/)) {
          // New feature
          if (currentFeature) {
            features.push(currentFeature);
          }
          
          const match = line.match(/^\s+(\S+)\s+(.*)/);
          if (match) {
            const type = match[1];
            const location = match[2];
            currentFeature = { type, location, qualifiers: {} };
          }
        } else if (line.match(/^\s+\/\S/) && currentFeature) {
          // Qualifier
          const match = line.trim().match(/\/(\S+)=(?:"([^"]*)"|([\S]+))/);
          if (match) {
            const key = match[1];
            const value = match[2] || match[3];
            currentFeature.qualifiers[key] = value;
          }
        }
      }
      
      if (currentFeature) {
        features.push(currentFeature);
      }
    }
  }
  
  // Count feature types
  const typeCounts = {};
  for (const feature of features) {
    typeCounts[feature.type] = (typeCounts[feature.type] || 0) + 1;
  }
  
  // Convert to array and sort by count
  const featureTypes = Object.entries(typeCounts)
    .map(([type, count]) => ({ type, count }))
    .sort((a, b) => b.count - a.count);
  
  // Get top examples of each feature type
  const topFeatures = [];
  const seenTypes = new Set();
  
  for (const feature of features) {
    if (!seenTypes.has(feature.type)) {
      topFeatures.push(feature);
      seenTypes.add(feature.type);
    }
    
    if (topFeatures.length >= 10) break;
  }
  
  return {
    metadata,
    featureCount: features.length,
    featureTypes,
    topFeatures,
    note: 'Feature locations may be simplified for display'
  };
};

/**
 * Analyze protein structure from PDB
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeStructure = (content, format) => {
  if (format !== 'pdb') {
    throw new Error('Structure analysis is only supported for PDB format');
  }
  
  const lines = content.split('\n');
  const atoms = [];
  const metadata = {};
  const chains = new Set();
  const residueTypes = new Map();
  
  // Extract header info
  const headerMatch = content.match(/HEADER\s+(.*)/);
  if (headerMatch) {
    metadata.header = headerMatch[1].trim();
  }
  
  const titleMatch = content.match(/TITLE\s+(.*?)(?=\nREMARK|\nSOURCE|\nAUTHOR|\nATOM)/s);
  if (titleMatch) {
    metadata.title = titleMatch[1].replace(/\n\s+/g, ' ').trim();
  }
  
  // Process atoms
  for (const line of lines) {
    if (line.startsWith('ATOM')) {
      const chainId = line.substring(21, 22).trim();
      const residueType = line.substring(17, 20).trim();
      
      // Keep track of chains
      chains.add(chainId);
      
      // Count residue types
      residueTypes.set(residueType, (residueTypes.get(residueType) || 0) + 1);
      
      // Store atom info
      atoms.push({
        atomNum: line.substring(6, 11).trim(),
        atomName: line.substring(12, 16).trim(),
        residueName: residueType,
        chainId: chainId,
        residueNum: line.substring(22, 26).trim(),
        x: parseFloat(line.substring(30, 38).trim()),
        y: parseFloat(line.substring(38, 46).trim()),
        z: parseFloat(line.substring(46, 54).trim())
      });
    }
  }
  
  // Extract secondary structure
  let helixCount = 0;
  let sheetCount = 0;
  
  for (const line of lines) {
    if (line.startsWith('HELIX')) {
      helixCount++;
    } else if (line.startsWith('SHEET')) {
      sheetCount++;
    }
  }
  
  // Sort residue types by count
  const residueComposition = Array.from(residueTypes.entries())
    .map(([residue, count]) => ({ residue, count }))
    .sort((a, b) => b.count - a.count);
  
  // Chain information
  const chainInfo = Array.from(chains).map(chainId => {
    // Count residues in this chain
    const chainResidues = new Set();
    for (const atom of atoms) {
      if (atom.chainId === chainId && atom.atomName === 'CA') {
        chainResidues.add(`${atom.residueName}_${atom.residueNum}`);
      }
    }
    
    return {
      id: chainId,
      residueCount: chainResidues.size,
      description: `Chain ${chainId} (${chainResidues.size} residues)`
    };
  }).sort((a, b) => b.residueCount - a.residueCount);
  
  return {
    summary: {
      title: metadata.title || metadata.header || 'Unknown structure',
      atomCount: atoms.length,
      hetatmCount: 0, // We're only counting ATOM records here
      chainCount: chains.size,
      chains: Array.from(chains),
      residueTypeCount: residueTypes.size,
      secondaryStructure: {
        helixCount,
        sheetCount
      }
    },
    residueComposition,
    chains: chainInfo,
    note: 'Basic structure analysis without calculating complex structural properties'
  };
};

/**
 * Analyze residue distribution in PDB
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeResidueDistribution = (content, format) => {
  if (format !== 'pdb') {
    throw new Error('Residue distribution analysis is only supported for PDB format');
  }
  
  const lines = content.split('\n');
  const residues = new Map();
  const chains = new Set();
  
  // Process atoms
  for (const line of lines) {
    if (line.startsWith('ATOM')) {
      const chainId = line.substring(21, 22).trim();
      const residueType = line.substring(17, 20).trim();
      const residueNum = line.substring(22, 26).trim();
      const atomName = line.substring(12, 16).trim();
      
      // Only count CA atoms to avoid counting residues multiple times
      if (atomName === 'CA') {
        const residueKey = `${chainId}_${residueNum}`;
        residues.set(residueKey, residueType);
        chains.add(chainId);
      }
    }
  }
  
  // Count residue types
  const residueCounts = {};
  for (const residueType of residues.values()) {
    residueCounts[residueType] = (residueCounts[residueType] || 0) + 1;
  }
  
  // Calculate percentages
  const totalResidues = residues.size;
  const residueStats = Object.entries(residueCounts)
    .map(([residue, count]) => ({
      residue,
      count,
      percentage: ((count / totalResidues) * 100).toFixed(2) + '%'
    }))
    .sort((a, b) => b.count - a.count);
  
  // Analyze chains
  const chainStats = Array.from(chains).map(chainId => {
    // Count residues in this chain
    let count = 0;
    for (const [key, _] of residues.entries()) {
      if (key.startsWith(chainId + '_')) {
        count++;
      }
    }
    
    return {
      chain: chainId,
      residueCount: count,
      percentage: ((count / totalResidues) * 100).toFixed(2) + '%'
    };
  }).sort((a, b) => b.residueCount - a.residueCount);
  
  // Group residues by type
  const residueGroups = {
    hydrophobic: ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'],
    polar: ['GLY', 'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'],
    charged: ['LYS', 'ARG', 'HIS', 'ASP', 'GLU']
  };
  
  const groupStats = {};
  for (const [group, residueList] of Object.entries(residueGroups)) {
    let count = 0;
    for (const residue of residueList) {
      count += residueCounts[residue] || 0;
    }
    
    groupStats[group] = {
      count,
      percentage: ((count / totalResidues) * 100).toFixed(2) + '%'
    };
  }
  
  return {
    totalResidues,
    residueStats,
    chainStats,
    groupStats,
    note: 'Residue counts are based on CA atoms to avoid counting residues multiple times'
  };
};

/**
 * Analyze conservation in multiple sequence alignments
 * @param {string} content - Content to analyze
 * @param {string} format - Format of the content
 * @returns {object} Analysis results
 */
const analyzeConservation = (content, format) => {
  if (format !== 'clustal' && format !== 'stockholm') {
    throw new Error('Conservation analysis is only supported for CLUSTAL and Stockholm formats');
  }
  
  const sequences = {};
  let conservation = '';
  
  if (format === 'clustal') {
    // Parse CLUSTAL format
    const lines = content.split('\n');
    
    // Find the line with CLUSTAL header
    const headerIndex = lines.findIndex(line => line.includes('CLUSTAL'));
    
    if (headerIndex !== -1) {
      // Process sequence and conservation lines
      for (let i = headerIndex + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        
        if (!line) continue;
        
        if (line.startsWith(' ')) {
          // Conservation line
          conservation += line;
        } else {
          // Sequence line
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
    }
  } else if (format === 'stockholm') {
    // Parse Stockholm format
    const lines = content.split('\n');
    
    for (const line of lines) {
      if (line.startsWith('#=GC') && line.includes('cons')) {
        // Conservation line
        const parts = line.split(/\s+/);
        if (parts.length >= 3) {
          conservation += parts[parts.length - 1];
        }
      } else if (!line.startsWith('#') && !line.startsWith('//') && line.trim()) {
        // Sequence line
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
    }
    
    // If no conservation line was found, generate one
    if (!conservation && Object.values(sequences).length > 0) {
      const firstSeq = Object.values(sequences)[0];
      conservation = generateConservation(Object.values(sequences));
    }
  }
  
  // Clean up conservation string
  conservation = conservation.replace(/\s+/g, '');
  
  // Calculate conservation statistics
  const conservationStats = analyzeConservationString(conservation);
  
  // Find conserved regions
  const conservedRegions = findConservedRegions(conservation);
  
  // Prepare visualization
  const visualization = createConservationVisualization(sequences, conservation);
  
  return {
    sequenceCount: Object.keys(sequences).length,
    alignmentLength: Object.values(sequences)[0]?.length || 0,
    overallConservation: conservationStats.percentage,
    conservationScore: conservationStats,
    conservedRegions,
    visualization
  };
};

/**
 * Generate a conservation string for aligned sequences
 * @param {Array} sequences - Array of sequence strings
 * @returns {string} Conservation string
 */
const generateConservation = (sequences) => {
  if (!sequences.length || !sequences[0].length) return '';
  
  let conservation = '';
  
  for (let i = 0; i < sequences[0].length; i++) {
    // Get all characters at this position
    const column = sequences.map(seq => i < seq.length ? seq[i] : '-');
    
    // Check if all characters are the same
    const allMatch = column.every(char => char === column[0] && char !== '-');
    
    // Mark conservation
    if (allMatch) {
      conservation += '*';
    } else {
      // Check if the column is somewhat conserved (e.g., similar amino acids)
      const uniqueChars = new Set(column.filter(char => char !== '-'));
      
      if (uniqueChars.size <= 2) {
        conservation += ':';
      } else if (uniqueChars.size <= 3) {
        conservation += '.';
      } else {
        conservation += ' ';
      }
    }
  }
  
  return conservation;
};

/**
 * Analyze a conservation string
 * @param {string} conservation - Conservation string
 * @returns {object} Conservation statistics
 */
const analyzeConservationString = (conservation) => {
  if (!conservation) return { percentage: 0 };
  
  const counts = {
    '*': 0, // Fully conserved
    ':': 0, // Strongly similar
    '.': 0, // Weakly similar
    ' ': 0  // Not conserved
  };
  
  // Count symbols
  for (const char of conservation) {
    counts[char] = (counts[char] || 0) + 1;
  }
  
  // Calculate percentages
  const total = conservation.length;
  const fullyConserved = counts['*'] || 0;
  const stronglySimilar = counts[':'] || 0;
  const weaklySimilar = counts['.'] || 0;
  
  // Overall conservation score (weighting: * = 1, : = 0.5, . = 0.25)
  const conservationScore = (
    (fullyConserved * 1.0) + 
    (stronglySimilar * 0.5) + 
    (weaklySimilar * 0.25)
  ) / total;
  
  return {
    total,
    fullyConserved,
    stronglySimilar,
    weaklySimilar,
    notConserved: total - fullyConserved - stronglySimilar - weaklySimilar,
    percentage: Math.round(conservationScore * 100)
  };
};

/**
 * Find conserved regions in the alignment
 * @param {string} conservation - Conservation string
 * @returns {Array} Array of conserved regions
 */
const findConservedRegions = (conservation) => {
  const regions = [];
  let start = -1;
  
  // Find runs of conserved positions (*)
  for (let i = 0; i < conservation.length; i++) {
    if (conservation[i] === '*') {
      if (start === -1) {
        start = i;
      }
    } else {
      if (start !== -1) {
        // Found the end of a conserved region
        const length = i - start;
        
        if (length >= 3) { // Only record regions of at least 3 conserved positions
          regions.push({
            start: start + 1, // 1-based position
            end: i,
            length,
            conservation: '100%'
          });
        }
        
        start = -1;
      }
    }
  }
  
  // Check if we ended with a conserved region
  if (start !== -1) {
    const length = conservation.length - start;
    
    if (length >= 3) {
      regions.push({
        start: start + 1,
        end: conservation.length,
        length,
        conservation: '100%'
      });
    }
  }
  
  // Also find regions with strong similarity
  start = -1;
  let strongCount = 0;
  
  for (let i = 0; i < conservation.length; i++) {
    if (conservation[i] === '*' || conservation[i] === ':') {
      if (start === -1) {
        start = i;
      }
      if (conservation[i] === '*') strongCount++;
    } else {
      if (start !== -1) {
        // Found the end of a similar region
        const length = i - start;
        
        if (length >= 5 && start !== -1) { // Only record regions of at least 5 positions
          const conservationPercent = Math.round((strongCount / length) * 100);
          
          // Only add if conservation is at least 70%
          if (conservationPercent >= 70) {
            regions.push({
              start: start + 1,
              end: i,
              length,
              conservation: `${conservationPercent}%`
            });
          }
        }
        
        start = -1;
        strongCount = 0;
      }
    }
  }
  
  // Check if we ended with a similar region
  if (start !== -1) {
    const length = conservation.length - start;
    
    if (length >= 5) {
      const conservationPercent = Math.round((strongCount / length) * 100);
      
      if (conservationPercent >= 70) {
        regions.push({
          start: start + 1,
          end: conservation.length,
          length,
          conservation: `${conservationPercent}%`
        });
      }
    }
  }
  
  // Sort by conservation and length
  return regions.sort((a, b) => {
    const aPercent = parseInt(a.conservation);
    const bPercent = parseInt(b.conservation);
    
    if (aPercent !== bPercent) {
      return bPercent - aPercent;
    }
    
    return b.length - a.length;
  });
};

/**
 * Create a simple visualization of the conservation
 * @param {object} sequences - Object with sequence names and strings
 * @param {string} conservation - Conservation string
 * @returns {Array} Array of lines for visualization
 */
const createConservationVisualization = (sequences, conservation) => {
  const visualization = [];
  const seqArray = Object.entries(sequences);
  const lineLength = 60;
  
  // Get first sequence length
  const firstSeqLength = seqArray[0]?.[1].length || 0;
  
  // Create position header
  for (let i = 0; i < firstSeqLength; i += lineLength) {
    const end = Math.min(i + lineLength, firstSeqLength);
    
    // Position numbers
    let posLine = '         ';
    for (let j = i + 1; j <= end; j += 10) {
      posLine += j.toString().padEnd(10);
    }
    visualization.push(posLine);
    
    // Position markers
    let markerLine = '         ';
    for (let j = i; j < end; j++) {
      markerLine += ((j + 1) % 10 === 0) ? '|' : '.';
    }
    visualization.push(markerLine);
    
    // Sequences
    for (const [name, sequence] of seqArray) {
      const displayName = name.length > 8 ? name.substring(0, 8) : name.padEnd(8);
      const segment = sequence.substring(i, end);
      visualization.push(`${displayName} ${segment}`);
    }
    
    // Conservation line
    if (conservation) {
      const consSegment = conservation.substring(i, end);
      visualization.push(`Conserv.  ${consSegment}`);
    }
    
    // Add a blank line between blocks
    visualization.push('');
  }
  
  return visualization;
};
