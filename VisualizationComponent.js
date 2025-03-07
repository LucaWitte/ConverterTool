// VisualizationComponent.js - Provides visualization for different biological data formats

import React, { useState, useEffect, useRef, useMemo } from 'react';
import { Layers, BarChart2, Activity, FileText, ChevronDown, ZoomIn, ZoomOut, RefreshCw, Download, Settings } from 'lucide-react';
import * as d3 from 'd3';
import * as THREE from 'three';
import { formatSequenceForFasta } from './FormatConverters';
import { detectSequenceType } from './FormatDetection';


/**
 * Main Visualization Component that renders different visualizations based on format
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
export const VisualizationComponent = ({ 
  format, 
  content, 
  data, 
  width = '100%', 
  height = 500, 
  onClose 
}) => {
  const [activeTab, setActiveTab] = useState('default');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [settings, setSettings] = useState({
    colorScheme: 'default',
    showLabels: true,
    renderQuality: 'medium',
    fontSize: 12,
    showRuler: true,
    showLegend: true,
    lineHeight: 18
  });
  const [showSettings, setShowSettings] = useState(false);
  
  // Determine available visualizations based on format
  const getAvailableVisualizations = (format) => {
    switch (format) {
      case 'fasta':
        return [
          { id: 'sequence', label: 'Sequence View', icon: 'Layers' },
          { id: 'composition', label: 'Composition Analysis', icon: 'BarChart2' },
          { id: 'hydrophobicity', label: 'Hydrophobicity Plot', icon: 'Activity' }
        ];
      case 'fastq':
        return [
          { id: 'sequence', label: 'Sequence View', icon: 'Layers' },
          { id: 'quality', label: 'Quality Scores', icon: 'Activity' }
        ];
      case 'genbank':
      case 'embl':
        return [
          { id: 'sequence', label: 'Sequence View', icon: 'Layers' },
          { id: 'features', label: 'Feature Map', icon: 'Layers' },
          { id: 'circular', label: 'Circular View', icon: 'RefreshCw' }
        ];
      case 'pdb':
      case 'mmcif':
        return [
          { id: '3dstructure', label: '3D Structure', icon: 'Cube' },
          { id: 'sequence', label: 'Sequence View', icon: 'Layers' },
          { id: 'contacts', label: 'Contact Map', icon: 'Grid' }
        ];
      case 'clustal':
      case 'stockholm':
        return [
          { id: 'alignment', label: 'Alignment View', icon: 'Layers' },
          { id: 'conservation', label: 'Conservation', icon: 'BarChart2' }
        ];
      case 'newick':
      case 'nexus':
        return [
          { id: 'tree', label: 'Phylogenetic Tree', icon: 'GitBranch' }
        ];
      case 'bed':
      case 'gff':
      case 'vcf':
        return [
          { id: 'features', label: 'Feature Map', icon: 'Layers' }
        ];
      default:
        return [
          { id: 'text', label: 'Text View', icon: 'FileText' }
        ];
    }
  };
  
  const visualizations = getAvailableVisualizations(format);
  
  // Set first visualization as default
  useEffect(() => {
    if (visualizations.length > 0 && !activeTab) {
      setActiveTab(visualizations[0].id);
    }
  }, [format, visualizations, activeTab]);
  
  // Render the appropriate visualization component
  const renderVisualization = () => {
    setLoading(true);
    setError(null);
    
    try {
      switch (activeTab) {
        case 'sequence':
          return <SequenceViewer content={content} data={data} format={format} settings={settings} />;
        case 'composition':
          return <CompositionAnalysis content={content} data={data} format={format} settings={settings} />;
        case 'hydrophobicity':
          return <HydrophobicityPlot content={content} data={data} settings={settings} />;
        case 'quality':
          return <QualityScoreViewer content={content} data={data} settings={settings} />;
        case 'features':
          return <FeatureMapViewer content={content} data={data} format={format} settings={settings} />;
        case 'circular':
          return <CircularViewer content={content} data={data} settings={settings} />;
        case '3dstructure':
          return <StructureViewer content={content} data={data} settings={settings} />;
        case 'contacts':
          return <ContactMapViewer content={content} data={data} settings={settings} />;
        case 'alignment':
          return <AlignmentViewer content={content} data={data} settings={settings} />;
        case 'conservation':
          return <ConservationViewer content={content} data={data} settings={settings} />;
        case 'tree':
          return <PhylogeneticTreeViewer content={content} data={data} settings={settings} />;
        case 'text':
        default:
          return <TextViewer content={content} settings={settings} />;
      }
    } catch (e) {
      setError(`Error rendering visualization: ${e.message}`);
      return <ErrorDisplay error={e.message} />;
    } finally {
      setLoading(false);
    }
  };
  
  return (
    <div className="bg-white border rounded-lg shadow-md overflow-hidden" style={{ width, height }}>
      {/* Visualization tabs */}
      <div className="border-b bg-gray-50 p-2 flex justify-between items-center">
        <div className="flex overflow-x-auto hide-scrollbar">
          {visualizations.map(viz => (
            <button
              key={viz.id}
              className={`px-3 py-1 flex items-center ${activeTab === viz.id ? 'bg-blue-600 text-white rounded' : 'text-gray-600 hover:bg-gray-100 rounded'}`}
              onClick={() => setActiveTab(viz.id)}
            >
              {viz.icon === 'Layers' && <Layers className="w-4 h-4 mr-1" />}
              {viz.icon === 'BarChart2' && <BarChart2 className="w-4 h-4 mr-1" />}
              {viz.icon === 'Activity' && <Activity className="w-4 h-4 mr-1" />}
              {viz.icon === 'FileText' && <FileText className="w-4 h-4 mr-1" />}
              {viz.icon === 'RefreshCw' && <RefreshCw className="w-4 h-4 mr-1" />}
              <span className="text-sm whitespace-nowrap">{viz.label}</span>
            </button>
          ))}
        </div>
        
        <div className="flex items-center">
          <button
            className="p-1 text-gray-500 hover:bg-gray-100 rounded mr-1"
            onClick={() => setShowSettings(!showSettings)}
            title="Visualization Settings"
          >
            <Settings className="w-4 h-4" />
          </button>
          {onClose && (
            <button
              className="p-1 text-gray-500 hover:bg-gray-100 rounded"
              onClick={onClose}
              title="Close Visualization"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          )}
        </div>
      </div>
      
      {/* Settings panel */}
      {showSettings && (
        <div className="border-b p-3 bg-gray-50">
          <h3 className="text-sm font-medium mb-2">Visualization Settings</h3>
          <div className="grid grid-cols-2 gap-3">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Color Scheme</label>
              <select
                className="w-full text-sm p-1 border rounded"
                value={settings.colorScheme}
                onChange={(e) => setSettings({...settings, colorScheme: e.target.value})}
              >
                <option value="default">Default</option>
                <option value="hydrophobicity">Hydrophobicity</option>
                <option value="charge">Charge</option>
                <option value="nucleotide">Nucleotide</option>
                <option value="taylor">Taylor</option>
                <option value="clustal">Clustal</option>
              </select>
            </div>
            
            <div>
              <label className="block text-xs text-gray-600 mb-1">Font Size</label>
              <select
                className="w-full text-sm p-1 border rounded"
                value={settings.fontSize}
                onChange={(e) => setSettings({...settings, fontSize: parseInt(e.target.value)})}
              >
                <option value="10">Small</option>
                <option value="12">Medium</option>
                <option value="14">Large</option>
                <option value="16">X-Large</option>
              </select>
            </div>
            
            <div className="flex items-center">
              <input
                type="checkbox"
                id="showLabels"
                checked={settings.showLabels}
                onChange={(e) => setSettings({...settings, showLabels: e.target.checked})}
                className="mr-2"
              />
              <label htmlFor="showLabels" className="text-xs text-gray-600">Show Labels</label>
            </div>
            
            <div className="flex items-center">
              <input
                type="checkbox"
                id="showRuler"
                checked={settings.showRuler}
                onChange={(e) => setSettings({...settings, showRuler: e.target.checked})}
                className="mr-2"
              />
              <label htmlFor="showRuler" className="text-xs text-gray-600">Show Ruler</label>
            </div>
            
            <div className="flex items-center">
              <input
                type="checkbox"
                id="showLegend"
                checked={settings.showLegend}
                onChange={(e) => setSettings({...settings, showLegend: e.target.checked})}
                className="mr-2"
              />
              <label htmlFor="showLegend" className="text-xs text-gray-600">Show Legend</label>
            </div>
            
            <div>
              <label className="block text-xs text-gray-600 mb-1">Render Quality</label>
              <select
                className="w-full text-sm p-1 border rounded"
                value={settings.renderQuality}
                onChange={(e) => setSettings({...settings, renderQuality: e.target.value})}
              >
                <option value="low">Low (Faster)</option>
                <option value="medium">Medium</option>
                <option value="high">High (Slower)</option>
              </select>
            </div>
          </div>
        </div>
      )}
      
      {/* Visualization content */}
      <div className="flex-1 overflow-hidden relative" style={{ height: showSettings ? 'calc(100% - 98px)' : 'calc(100% - 42px)' }}>
        {loading ? (
          <div className="flex items-center justify-center h-full">
            <RefreshCw className="w-6 h-6 text-blue-600 animate-spin" />
            <span className="ml-2 text-gray-600">Loading visualization...</span>
          </div>
        ) : error ? (
          <ErrorDisplay error={error} />
        ) : (
          renderVisualization()
        )}
      </div>
    </div>
  );
};

/**
 * Error display component
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const ErrorDisplay = ({ error }) => (
  <div className="flex flex-col items-center justify-center h-full p-4 text-center">
    <svg className="w-12 h-12 text-red-500 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4m0 4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
    </svg>
    <h3 className="text-lg font-medium text-red-600">Visualization Error</h3>
    <p className="text-gray-600 mt-1">{error}</p>
  </div>
);

/**
 * Sequence viewer component for nucleotide/protein sequences
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const SequenceViewer = ({ content, data, format, settings }) => {
  const containerRef = useRef(null);
  const [sequences, setSequences] = useState([]);
  const [zoomLevel, setZoomLevel] = useState(1);
  const [viewportStart, setViewportStart] = useState(0);
  const [viewportEnd, setViewportEnd] = useState(0);
  const [sequenceType, setSequenceType] = useState('unknown');
  
  useEffect(() => {
    // Extract sequences based on format
    const extractedSequences = extractSequences(content, format);
    setSequences(extractedSequences);
    
    // Set viewport for the sequence
    if (extractedSequences.length > 0) {
      const firstSeq = extractedSequences[0];
      setViewportEnd(Math.min(200, firstSeq.sequence.length));
      
      // Detect sequence type
      const detectedType = detectSequenceType(firstSeq.sequence);
      setSequenceType(detectedType);
    }
  }, [content, format]);
  
  // Extract sequences from different formats
  const extractSequences = (content, format) => {
    if (!content) return [];
    
    switch (format) {
      case 'fasta':
        return extractFastaSequences(content);
      case 'genbank':
      case 'embl':
        return extractAnnotatedSequences(content, format);
      case 'pdb':
      case 'mmcif':
        return extractProteinSequences(content, format);
      case 'clustal':
      case 'stockholm':
        return extractAlignmentSequences(content, format);
      default:
        return [];
    }
  };
  
  // Extract FASTA sequences
  const extractFastaSequences = (content) => {
    const sequences = [];
    const lines = content.split('\n');
    
    let currentSeq = null;
    let currentHeader = "";
    
    for (const line of lines) {
      if (line.startsWith('>')) {
        if (currentSeq) {
          sequences.push({
            header: currentHeader,
            sequence: currentSeq,
            length: currentSeq.length
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
        sequence: currentSeq,
        length: currentSeq.length
      });
    }
    
    return sequences;
  };
  
  // Extract sequences from GenBank/EMBL
  const extractAnnotatedSequences = (content, format) => {
    // This would be a more complex parser for GenBank/EMBL
    // For now, simplified to extract just the sequence
    const sequences = [];
    
    if (format === 'genbank') {
      const originMatch = content.match(/ORIGIN([\s\S]*?)\/\//);
      if (originMatch) {
        let sequence = originMatch[1].replace(/\d+|\s+/g, '');
        
        const locusMatch = content.match(/LOCUS\s+(\S+)/);
        const header = locusMatch ? locusMatch[1] : 'Unknown';
        
        sequences.push({
          header,
          sequence,
          length: sequence.length
        });
      }
    } else if (format === 'embl') {
      const sequenceMatch = content.match(/SQ.*?Sequence([\s\S]*?)\/\//);
      if (sequenceMatch) {
        let sequence = sequenceMatch[1].replace(/\d+|\s+/g, '');
        
        const idMatch = content.match(/ID\s+(\S+)/);
        const header = idMatch ? idMatch[1] : 'Unknown';
        
        sequences.push({
          header,
          sequence,
          length: sequence.length
        });
      }
    }
    
    return sequences;
  };
  
  // Extract protein sequences from PDB/mmCIF
  const extractProteinSequences = (content, format) => {
    // Simplified extraction for demo
    const sequences = [];
    
    if (format === 'pdb') {
      // Extract chains from PDB
      const chainMap = new Map();
      const lines = content.split('\n');
      
      for (const line of lines) {
        if (line.startsWith('ATOM') && line.length >= 26) {
          const chainId = line.substring(21, 22).trim();
          const residue = line.substring(17, 20).trim();
          const residueNum = parseInt(line.substring(22, 26).trim());
          
          if (!chainMap.has(chainId)) {
            chainMap.set(chainId, new Map());
          }
          
          chainMap.get(chainId).set(residueNum, residue);
        }
      }
      
      // Convert 3-letter to 1-letter amino acids
      const aa3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
      };
      
      // Create sequences for each chain
      for (const [chainId, residueMap] of chainMap.entries()) {
        const sortedResidues = [...residueMap.entries()].sort((a, b) => a[0] - b[0]);
        
        let sequence = '';
        for (const [_, residue] of sortedResidues) {
          sequence += aa3to1[residue] || 'X';
        }
        
        sequences.push({
          header: `Chain ${chainId}`,
          sequence,
          length: sequence.length
        });
      }
    }
    
    return sequences;
  };
  
  // Extract sequences from alignments
  const extractAlignmentSequences = (content, format) => {
    const sequences = [];
    
    if (format === 'clustal') {
      const lines = content.split('\n');
      
      // Find the header line
      const headerLineIndex = lines.findIndex(line => line.includes('CLUSTAL'));
      
      if (headerLineIndex !== -1) {
        const sequenceMap = {};
        
        // Parse sequences
        for (let i = headerLineIndex + 1; i < lines.length; i++) {
          const line = lines[i].trim();
          if (!line || line.startsWith(' ')) continue;
          
          const parts = line.split(/\s+/);
          if (parts.length >= 2) {
            const seqName = parts[0];
            const seqSegment = parts[1];
            
            if (!sequenceMap[seqName]) {
              sequenceMap[seqName] = '';
            }
            
            sequenceMap[seqName] += seqSegment;
          }
        }
        
        // Convert to array
        for (const [name, sequence] of Object.entries(sequenceMap)) {
          sequences.push({
            header: name,
            sequence,
            length: sequence.length
          });
        }
      }
    }
    
    return sequences;
  };
  
  // Get color for residue based on property
  const getResidueColor = (residue, colorScheme, seqType) => {
    const residue_upper = residue.toUpperCase();
    
    // Default color
    if (!residue_upper || residue_upper === '-' || residue_upper === '.') {
      return '#cccccc';
    }
    
    const scheme = colorScheme || 'default';
    
    // Different color schemes
    const schemes = {
      default: {
        A: '#77dd88', C: '#99ee66', D: '#55bb33', E: '#55bb33',
        F: '#9999ff', G: '#77dd88', H: '#7070ff', I: '#66bbff',
        K: '#ffcc77', L: '#66bbff', M: '#66bbff', N: '#55bb33',
        P: '#eeaaaa', Q: '#55bb33', R: '#ffcc77', S: '#ff4455',
        T: '#ff4455', V: '#66bbff', W: '#9999ff', Y: '#9999ff'
      },
      hydrophobicity: {
        A: '#ad0052', C: '#0000ff', D: '#0000ff', E: '#0000ff',
        F: '#ad0052', G: '#0000ff', H: '#0000ff', I: '#ad0052',
        K: '#0000ff', L: '#ad0052', M: '#ad0052', N: '#0000ff',
        P: '#0000ff', Q: '#0000ff', R: '#0000ff', S: '#0000ff',
        T: '#0000ff', V: '#ad0052', W: '#ad0052', Y: '#0000ff'
      },
      charge: {
        A: '#000000', C: '#000000', D: '#ff0000', E: '#ff0000',
        F: '#000000', G: '#000000', H: '#0000ff', I: '#000000',
        K: '#0000ff', L: '#000000', M: '#000000', N: '#000000',
        P: '#000000', Q: '#000000', R: '#0000ff', S: '#000000',
        T: '#000000', V: '#000000', W: '#000000', Y: '#000000'
      },
      nucleotide: {
        A: '#a0a0ff', T: '#a0ffa0', G: '#ff7070', C: '#ffff70',
        U: '#a0ffa0', N: '#ffffff'
      },
      clustal: {
        // Clustal color scheme
        A: '#80a0f0', R: '#f01505', N: '#00ff00', D: '#c048c0',
        C: '#f08080', Q: '#00ff00', E: '#c048c0', G: '#f09048',
        H: '#15a4a4', I: '#80a0f0', L: '#80a0f0', K: '#f01505',
        M: '#80a0f0', F: '#80a0f0', P: '#ffff00', S: '#00ff00',
        T: '#00ff00', W: '#80a0f0', Y: '#15a4a4', V: '#80a0f0',
        B: '#fff', Z: '#fff', X: '#fff', '-': '#fff'
      }
    };
    
    if (seqType === 'dna' || seqType === 'rna') {
      return schemes.nucleotide[residue_upper] || '#cccccc';
    } else {
      return schemes[scheme][residue_upper] || '#cccccc';
    }
  };
  
  // Render
  const renderSequences = () => {
    if (sequences.length === 0) {
      return (
        <div className="flex items-center justify-center h-full">
          <p className="text-gray-500">No sequence data available</p>
        </div>
      );
    }
    
    const lineHeight = settings.lineHeight || 18;
    const lineWidth = 60; // Characters per line
    const fontSize = settings.fontSize || 12;
    
    return (
      <div className="p-4">
        {sequences.map((seq, index) => (
          <div key={index} className="mb-6">
            <div className="font-medium mb-2 text-blue-600">{seq.header}</div>
            
            {/* Render sequence */}
            <div 
              className="font-mono whitespace-pre" 
              style={{ fontSize: `${fontSize}px`, lineHeight: `${lineHeight}px` }}
            >
              {Array.from({ length: Math.ceil(seq.sequence.length / lineWidth) }).map((_, lineIndex) => {
                const start = lineIndex * lineWidth;
                const end = Math.min(start + lineWidth, seq.sequence.length);
                const lineSeq = seq.sequence.substring(start, end);
                
                return (
                  <div key={lineIndex} className="flex mb-1">
                    {/* Position number */}
                    {settings.showRuler && (
                      <div className="text-gray-500 w-12 text-right pr-2 select-none">
                        {start + 1}
                      </div>
                    )}
                    
                    {/* Sequence with colored residues */}
                    <div className="flex">
                      {Array.from(lineSeq).map((residue, i) => (
                        <span 
                          key={i} 
                          className="inline-block"
                          style={{ 
                            backgroundColor: getResidueColor(residue, settings.colorScheme, sequenceType),
                            width: `${fontSize * 0.8}px`,
                            height: `${fontSize * 1.2}px`,
                            textAlign: 'center',
                            color: settings.colorScheme === 'hydrophobicity' || settings.colorScheme === 'charge' ? 'white' : 'black',
                            lineHeight: `${fontSize * 1.2}px`
                          }}
                          title={`Position ${start + i + 1}: ${residue}`}
                        >
                          {residue}
                        </span>
                      ))}
                    </div>
                  </div>
                );
              })}
            </div>
          </div>
        ))}
      </div>
    );
  };
  
  return (
    <div className="h-full overflow-auto" ref={containerRef}>
      {/* Controls */}
      <div className="sticky top-0 bg-white border-b flex items-center p-2 z-10">
        <div className="flex-1"></div>
        <div className="flex items-center space-x-1">
          <button
            className="p-1 text-gray-500 hover:bg-gray-100 rounded"
            onClick={() => setZoomLevel(Math.max(0.5, zoomLevel - 0.1))}
            title="Zoom Out"
          >
            <ZoomOut className="w-4 h-4" />
          </button>
          <div className="text-sm text-gray-600">{Math.round(zoomLevel * 100)}%</div>
          <button
            className="p-1 text-gray-500 hover:bg-gray-100 rounded"
            onClick={() => setZoomLevel(Math.min(3, zoomLevel + 0.1))}
            title="Zoom In"
          >
            <ZoomIn className="w-4 h-4" />
          </button>
        </div>
      </div>
      
      {/* Sequence display */}
      <div className="p-2" style={{ transform: `scale(${zoomLevel})`, transformOrigin: 'top left' }}>
        {renderSequences()}
      </div>
    </div>
  );
};

/**
 * Composition analysis component
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const CompositionAnalysis = ({ content, data, format, settings }) => {
  const chartRef = useRef(null);
  const [sequences, setSequences] = useState([]);
  const [composition, setComposition] = useState({});
  const [activeChart, setActiveChart] = useState('frequency');
  
  useEffect(() => {
    // Extract sequences
    if (format === 'fasta') {
      const extractedSequences = extractFastaSequences(content);
      setSequences(extractedSequences);
      
      // Calculate composition
      if (extractedSequences.length > 0) {
        const seqType = detectSequenceType(extractedSequences[0].sequence);
        calculateComposition(extractedSequences, seqType);
      }
    }
  }, [content, format]);
  
  // Extract FASTA sequences
  const extractFastaSequences = (content) => {
    const sequences = [];
    const lines = content.split('\n');
    
    let currentSeq = null;
    let currentHeader = "";
    
    for (const line of lines) {
      if (line.startsWith('>')) {
        if (currentSeq) {
          sequences.push({
            header: currentHeader,
            sequence: currentSeq,
            length: currentSeq.length
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
        sequence: currentSeq,
        length: currentSeq.length
      });
    }
    
    return sequences;
  };
  
  // Calculate composition statistics
  const calculateComposition = (sequences, seqType) => {
    const isNucleotide = seqType === 'dna' || seqType === 'rna';
    
    // Define the elements to count based on sequence type
    const elements = isNucleotide 
      ? ['A', 'T', 'G', 'C', 'U', 'N'] 
      : ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'];
    
    // Initialize counts
    const counts = {};
    elements.forEach(e => { counts[e] = 0; });
    
    let totalChars = 0;
    
    // Count occurrences
    for (const seq of sequences) {
      const upper = seq.sequence.toUpperCase();
      for (const char of upper) {
        if (char !== '-' && char !== ' ') {
          counts[char] = (counts[char] || 0) + 1;
          totalChars++;
        }
      }
    }
    
    // Calculate percentages
    const percentages = {};
    for (const element of elements) {
      if (counts[element]) {
        percentages[element] = (counts[element] / totalChars) * 100;
      }
    }
    
    // Codon usage (for nucleotide sequences)
    let codonUsage = {};
    if (isNucleotide) {
      codonUsage = calculateCodonUsage(sequences);
    }
    
    // Amino acid properties (for protein sequences)
    let aminoAcidProps = {};
    if (!isNucleotide) {
      aminoAcidProps = calculateAminoAcidProps(sequences);
    }
    
    setComposition({
      counts,
      percentages,
      totalChars,
      seqType,
      codonUsage,
      aminoAcidProps
    });
  };
  
  // Calculate codon usage for nucleotide sequences
  const calculateCodonUsage = (sequences) => {
    const codonCounts = {};
    let totalCodons = 0;
    
    for (const seq of sequences) {
      const upper = seq.sequence.toUpperCase();
      for (let i = 0; i < upper.length - 2; i += 3) {
        const codon = upper.substring(i, i + 3);
        if (codon.length === 3 && !/[^ATGCUN]/.test(codon)) {
          codonCounts[codon] = (codonCounts[codon] || 0) + 1;
          totalCodons++;
        }
      }
    }
    
    // Calculate percentages
    const percentages = {};
    for (const [codon, count] of Object.entries(codonCounts)) {
      percentages[codon] = (count / totalCodons) * 100;
    }
    
    return { counts: codonCounts, percentages, totalCodons };
  };
  
  // Calculate amino acid properties for protein sequences
  const calculateAminoAcidProps = (sequences) => {
    // Define properties
    const properties = {
      hydrophobic: ['A', 'I', 'L', 'M', 'F', 'V', 'P', 'G'],
      polar: ['S', 'T', 'Y', 'C', 'N', 'Q'],
      charged: ['D', 'E', 'K', 'R', 'H'],
      positive: ['K', 'R', 'H'],
      negative: ['D', 'E']
    };
    
    // Count by property
    const propertyCounts = {};
    Object.keys(properties).forEach(prop => {
      propertyCounts[prop] = 0;
    });
    
    let totalAA = 0;
    
    // Count amino acids by property
    for (const seq of sequences) {
      const upper = seq.sequence.toUpperCase();
      for (const char of upper) {
        if (char !== '-' && char !== ' ' && char !== 'X') {
          totalAA++;
          
          // Add to appropriate property groups
          for (const [prop, aas] of Object.entries(properties)) {
            if (aas.includes(char)) {
              propertyCounts[prop]++;
            }
          }
        }
      }
    }
    
    // Calculate percentages
    const percentages = {};
    for (const [prop, count] of Object.entries(propertyCounts)) {
      percentages[prop] = (count / totalAA) * 100;
    }
    
    return { counts: propertyCounts, percentages, totalAA };
  };
  
  // Render frequency chart
  useEffect(() => {
    if (!chartRef.current || !composition.percentages) return;
    
    // Clear previous chart
    d3.select(chartRef.current).selectAll("*").remove();
    
    // Prepare data
    let chartData = [];
    
    if (activeChart === 'frequency') {
      chartData = Object.entries(composition.percentages)
        .filter(([_, value]) => value > 0)
        .map(([key, value]) => ({ 
          label: key, 
          value: value,
          count: composition.counts[key] || 0
        }))
        .sort((a, b) => b.value - a.value);
    } else if (activeChart === 'properties' && composition.seqType !== 'dna' && composition.seqType !== 'rna') {
      chartData = Object.entries(composition.aminoAcidProps.percentages)
        .map(([key, value]) => ({ 
          label: key, 
          value: value,
          count: composition.aminoAcidProps.counts[key] || 0
        }))
        .sort((a, b) => b.value - a.value);
    } else if (activeChart === 'codon' && (composition.seqType === 'dna' || composition.seqType === 'rna')) {
      // Top 20 codons
      chartData = Object.entries(composition.codonUsage.percentages)
        .filter(([_, value]) => value > 0)
        .map(([key, value]) => ({ 
          label: key, 
          value: value,
          count: composition.codonUsage.counts[key] || 0
        }))
        .sort((a, b) => b.value - a.value)
        .slice(0, 20);
    }
    
    // Set up dimensions
    const margin = { top: 30, right: 30, bottom: 70, left: 60 };
    const width = chartRef.current.clientWidth - margin.left - margin.right;
    const height = chartRef.current.clientHeight - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(chartRef.current)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // X axis
    const x = d3.scaleBand()
      .range([0, width])
      .domain(chartData.map(d => d.label))
      .padding(0.2);
    
    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(x))
      .selectAll("text")
        .attr("transform", "translate(-10,0)rotate(-45)")
        .style("text-anchor", "end");
    
    // Y axis
    const y = d3.scaleLinear()
      .domain([0, d3.max(chartData, d => d.value) * 1.1])
      .range([height, 0]);
    
    svg.append("g")
      .call(d3.axisLeft(y));
    
    // Add Y axis label
    svg.append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 0 - margin.left)
      .attr("x", 0 - (height / 2))
      .attr("dy", "1em")
      .style("text-anchor", "middle")
      .text("Percentage (%)");
    
    // Color scale
    const getBarColor = (label) => {
      if (composition.seqType === 'dna' || composition.seqType === 'rna') {
        const nucleotideColors = {
          'A': '#a0a0ff', 'T': '#a0ffa0', 'G': '#ff7070', 'C': '#ffff70',
          'U': '#a0ffa0', 'N': '#ffffff'
        };
        return nucleotideColors[label] || '#cccccc';
      } else {
        const aaColors = {
          'A': '#77dd88', 'C': '#99ee66', 'D': '#55bb33', 'E': '#55bb33',
          'F': '#9999ff', 'G': '#77dd88', 'H': '#7070ff', 'I': '#66bbff',
          'K': '#ffcc77', 'L': '#66bbff', 'M': '#66bbff', 'N': '#55bb33',
          'P': '#eeaaaa', 'Q': '#55bb33', 'R': '#ffcc77', 'S': '#ff4455',
          'T': '#ff4455', 'V': '#66bbff', 'W': '#9999ff', 'Y': '#9999ff',
          
          // Property colors
          'hydrophobic': '#66bbff',
          'polar': '#55bb33',
          'charged': '#ffcc77',
          'positive': '#ff7070',
          'negative': '#7070ff'
        };
        return aaColors[label] || '#cccccc';
      }
    };
    
    // Tooltip
    const tooltip = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0)
      .style("background-color", "white")
      .style("border", "solid")
      .style("border-width", "1px")
      .style("border-radius", "5px")
      .style("padding", "10px")
      .style("position", "absolute");
    
    // Add bars
    svg.selectAll("rect")
      .data(chartData)
      .join("rect")
        .attr("x", d => x(d.label))
        .attr("y", d => y(d.value))
        .attr("width", x.bandwidth())
        .attr("height", d => height - y(d.value))
        .attr("fill", d => getBarColor(d.label))
        .on("mouseover", function(d) {
          tooltip.transition()
            .duration(200)
            .style("opacity", .9);
          tooltip.html(`${d.label}: ${d.value.toFixed(2)}%<br>Count: ${d.count}`)
            .style("left", (event.pageX) + "px")
            .style("top", (event.pageY - 28) + "px");
        })
        .on("mouseout", function(d) {
          tooltip.transition()
            .duration(500)
            .style("opacity", 0);
        });
    
    // Chart title
    svg.append("text")
      .attr("x", width / 2)
      .attr("y", 0 - (margin.top / 2))
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .text(activeChart === 'frequency' ? 'Residue Frequency' : 
           activeChart === 'properties' ? 'Amino Acid Properties' : 'Codon Usage');
    
    // Return cleanup function
    return () => {
      d3.select("body").selectAll(".tooltip").remove();
    };
  }, [composition, activeChart]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Tab selection */}
      <div className="border-b flex">
        <button
          className={`px-4 py-2 ${activeChart === 'frequency' ? 'bg-blue-500 text-white' : 'hover:bg-gray-100'}`}
          onClick={() => setActiveChart('frequency')}
        >
          Residue Frequency
        </button>
        
        {composition.seqType === 'dna' || composition.seqType === 'rna' ? (
          <button
            className={`px-4 py-2 ${activeChart === 'codon' ? 'bg-blue-500 text-white' : 'hover:bg-gray-100'}`}
            onClick={() => setActiveChart('codon')}
          >
            Codon Usage
          </button>
        ) : (
          <button
            className={`px-4 py-2 ${activeChart === 'properties' ? 'bg-blue-500 text-white' : 'hover:bg-gray-100'}`}
            onClick={() => setActiveChart('properties')}
          >
            AA Properties
          </button>
        )}
      </div>
      
      {/* Chart */}
      <div className="flex-1 p-2">
        {composition.percentages ? (
          <div className="h-full" ref={chartRef}></div>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No composition data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Hydrophobicity plot for protein sequences
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const HydrophobicityPlot = ({ content, data, settings }) => {
  const chartRef = useRef(null);
  const [sequences, setSequences] = useState([]);
  const [activeSequence, setActiveSequence] = useState(0);
  const [windowSize, setWindowSize] = useState(9);
  
  useEffect(() => {
    // Extract sequences
    if (content) {
      const extractedSequences = extractFastaSequences(content);
      setSequences(extractedSequences);
    }
  }, [content]);
  
  // Extract FASTA sequences
  const extractFastaSequences = (content) => {
    const sequences = [];
    const lines = content.split('\n');
    
    let currentSeq = null;
    let currentHeader = "";
    
    for (const line of lines) {
      if (line.startsWith('>')) {
        if (currentSeq) {
          sequences.push({
            header: currentHeader,
            sequence: currentSeq,
            length: currentSeq.length
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
        sequence: currentSeq,
        length: currentSeq.length
      });
    }
    
    return sequences;
  };
  
  // Calculate hydrophobicity profile
  const calculateHydrophobicity = (sequence, windowSize) => {
    // Kyte & Doolittle hydrophobicity scale
    const scale = {
      'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
      'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
      'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
      'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    };
    
    const profile = [];
    
    for (let i = 0; i < sequence.length - windowSize + 1; i++) {
      let sum = 0;
      
      for (let j = 0; j < windowSize; j++) {
        const aa = sequence[i + j].toUpperCase();
        sum += scale[aa] || 0;
      }
      
      const avg = sum / windowSize;
      profile.push({ position: i + Math.floor(windowSize/2) + 1, value: avg });
    }
    
    return profile;
  };
  
  // Render hydrophobicity plot
  useEffect(() => {
    if (!chartRef.current || sequences.length === 0 || activeSequence >= sequences.length) return;
    
    const sequence = sequences[activeSequence].sequence;
    const profile = calculateHydrophobicity(sequence, windowSize);
    
    // Clear previous chart
    d3.select(chartRef.current).selectAll("*").remove();
    
    // Set up dimensions
    const margin = { top: 30, right: 30, bottom: 40, left: 60 };
    const width = chartRef.current.clientWidth - margin.left - margin.right;
    const height = chartRef.current.clientHeight - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(chartRef.current)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // X axis
    const x = d3.scaleLinear()
      .domain([1, sequence.length])
      .range([0, width]);
    
    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(x).ticks(10))
      .append("text")
        .attr("x", width / 2)
        .attr("y", margin.bottom - 5)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Residue Position");
    
    // Y axis
    const y = d3.scaleLinear()
      .domain([d3.min(profile, d => d.value) - 0.5, d3.max(profile, d => d.value) + 0.5])
      .range([height, 0]);
    
    svg.append("g")
      .call(d3.axisLeft(y))
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left + 15)
        .attr("x", -height / 2)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Hydrophobicity");
    
    // Add horizontal line at y=0
    svg.append("line")
      .attr("x1", 0)
      .attr("y1", y(0))
      .attr("x2", width)
      .attr("y2", y(0))
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .attr("stroke-dasharray", "4");
    
    // Add the line
    svg.append("path")
      .datum(profile)
      .attr("fill", "none")
      .attr("stroke", "#69b3a2")
      .attr("stroke-width", 2)
      .attr("d", d3.line()
        .x(d => x(d.position))
        .y(d => y(d.value))
      );
    
    // Add the area
    svg.append("path")
      .datum(profile)
      .attr("fill", "#69b3a2")
      .attr("fill-opacity", 0.3)
      .attr("stroke", "none")
      .attr("d", d3.area()
        .x(d => x(d.position))
        .y0(y(0))
        .y1(d => y(d.value))
      );
    
    // Add title
    svg.append("text")
      .attr("x", width / 2)
      .attr("y", -10)
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .text(`Hydrophobicity Plot (Window Size: ${windowSize})`);
    
    // Add tooltip
    const tooltip = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0)
      .style("background-color", "white")
      .style("border", "solid")
      .style("border-width", "1px")
      .style("border-radius", "5px")
      .style("padding", "5px")
      .style("position", "absolute");
    
    // Add dots and interaction
    svg.selectAll("circle")
      .data(profile)
      .join("circle")
        .attr("cx", d => x(d.position))
        .attr("cy", d => y(d.value))
        .attr("r", 0) // Start with radius 0
        .attr("fill", d => d.value > 0 ? "#ff7f0e" : "#1f77b4")
        .on("mouseover", function(d) {
          d3.select(this).transition()
            .duration(200)
            .attr("r", 5);
            
          tooltip.transition()
            .duration(200)
            .style("opacity", .9);
          tooltip.html(`Position: ${d.position}<br>Hydrophobicity: ${d.value.toFixed(2)}`)
            .style("left", (event.pageX + 10) + "px")
            .style("top", (event.pageY - 28) + "px");
        })
        .on("mouseout", function(d) {
          d3.select(this).transition()
            .duration(200)
            .attr("r", 0);
            
          tooltip.transition()
            .duration(500)
            .style("opacity", 0);
        });
    
    // Return cleanup function
    return () => {
      d3.select("body").selectAll(".tooltip").remove();
    };
  }, [sequences, activeSequence, windowSize]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div className="flex items-center space-x-2">
          <label className="text-sm text-gray-600">Sequence:</label>
          <select
            className="border p-1 rounded text-sm"
            value={activeSequence}
            onChange={(e) => setActiveSequence(parseInt(e.target.value))}
          >
            {sequences.map((seq, i) => (
              <option key={i} value={i}>
                {seq.header.length > 30 ? seq.header.substring(0, 30) + '...' : seq.header}
              </option>
            ))}
          </select>
        </div>
        
        <div className="flex items-center space-x-2">
          <label className="text-sm text-gray-600">Window Size:</label>
          <input
            type="range"
            min="3"
            max="21"
            step="2"
            value={windowSize}
            onChange={(e) => setWindowSize(parseInt(e.target.value))}
            className="w-24"
          />
          <span className="text-sm">{windowSize}</span>
        </div>
      </div>
      
      {/* Plot */}
      <div className="flex-1 p-2">
        {sequences.length > 0 ? (
          <div className="h-full" ref={chartRef}></div>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No sequence data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Quality score viewer for FASTQ files
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const QualityScoreViewer = ({ content, data, settings }) => {
  const chartRef = useRef(null);
  const [reads, setReads] = useState([]);
  const [stats, setStats] = useState(null);
  
  useEffect(() => {
    // Extract reads from FASTQ
    if (content) {
      const extractedReads = extractFastqReads(content);
      setReads(extractedReads);
      
      // Calculate quality statistics
      if (extractedReads.length > 0) {
        const qualityStats = calculateQualityStats(extractedReads);
        setStats(qualityStats);
      }
    }
  }, [content]);
  
  // Extract reads from FASTQ
  const extractFastqReads = (content) => {
    const reads = [];
    const lines = content.split('\n');
    
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
    
    return reads;
  };
  
  // Calculate quality statistics
  const calculateQualityStats = (reads) => {
    // Quality score by position
    const positionStats = {};
    
    // Get a representative sample (max 100 reads)
    const sampleSize = Math.min(reads.length, 100);
    const sample = reads.slice(0, sampleSize);
    
    // Find max read length
    const maxReadLength = Math.max(...sample.map(read => read.length));
    
    // Initialize position stats
    for (let i = 0; i < maxReadLength; i++) {
      positionStats[i] = {
        scores: [],
        mean: 0,
        median: 0,
        min: Infinity,
        max: -Infinity,
        q1: 0,
        q3: 0
      };
    }
    
    // Collect scores by position
    for (const read of sample) {
      for (let i = 0; i < read.quality.length; i++) {
        const qualScore = read.quality.charCodeAt(i) - 33; // Convert from ASCII to Phred
        positionStats[i].scores.push(qualScore);
      }
    }
    
    // Calculate statistics for each position
    for (let i = 0; i < maxReadLength; i++) {
      const scores = positionStats[i].scores.sort((a, b) => a - b);
      const count = scores.length;
      
      if (count > 0) {
        positionStats[i].min = scores[0];
        positionStats[i].max = scores[count - 1];
        positionStats[i].mean = scores.reduce((sum, score) => sum + score, 0) / count;
        positionStats[i].median = count % 2 === 0 ? 
          (scores[count/2 - 1] + scores[count/2]) / 2 : 
          scores[Math.floor(count/2)];
        positionStats[i].q1 = scores[Math.floor(count/4)];
        positionStats[i].q3 = scores[Math.floor(3*count/4)];
      }
    }
    
    // Convert to array format for d3
    const positionData = Object.entries(positionStats).map(([pos, stats]) => ({
      position: parseInt(pos) + 1,
      mean: stats.mean,
      median: stats.median,
      min: stats.min,
      max: stats.max,
      q1: stats.q1,
      q3: stats.q3
    }));
    
    // Overall stats
    const allScores = [];
    for (const read of sample) {
      for (let i = 0; i < read.quality.length; i++) {
        allScores.push(read.quality.charCodeAt(i) - 33);
      }
    }
    
    const sortedScores = allScores.sort((a, b) => a - b);
    
    const overallStats = {
      min: sortedScores[0],
      max: sortedScores[sortedScores.length - 1],
      mean: sortedScores.reduce((sum, score) => sum + score, 0) / sortedScores.length,
      median: sortedScores.length % 2 === 0 ?
        (sortedScores[sortedScores.length/2 - 1] + sortedScores[sortedScores.length/2]) / 2 :
        sortedScores[Math.floor(sortedScores.length/2)]
    };
    
    return {
      positionData,
      overallStats,
      readCount: reads.length,
      maxReadLength
    };
  };
  
  // Render quality plot
  useEffect(() => {
    if (!chartRef.current || !stats || !stats.positionData || stats.positionData.length === 0) return;
    
    // Clear previous chart
    d3.select(chartRef.current).selectAll("*").remove();
    
    // Set up dimensions
    const margin = { top: 30, right: 30, bottom: 40, left: 50 };
    const width = chartRef.current.clientWidth - margin.left - margin.right;
    const height = chartRef.current.clientHeight - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(chartRef.current)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // X axis
    const x = d3.scaleLinear()
      .domain([1, stats.maxReadLength])
      .range([0, width]);
    
    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(x).ticks(10))
      .append("text")
        .attr("x", width / 2)
        .attr("y", margin.bottom - 5)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Position in Read");
    
    // Y axis
    const y = d3.scaleLinear()
      .domain([0, 42]) // Phred quality score typically ranges from 0 to ~42
      .range([height, 0]);
    
    svg.append("g")
      .call(d3.axisLeft(y))
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left + 15)
        .attr("x", -height / 2)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Quality Score (Phred)");
    
    // Add quality zones
    svg.append("rect")
      .attr("x", 0)
      .attr("y", y(41))
      .attr("width", width)
      .attr("height", y(30) - y(41))
      .attr("fill", "#c3e6c3")
      .attr("opacity", 0.3);
    
    svg.append("rect")
      .attr("x", 0)
      .attr("y", y(30))
      .attr("width", width)
      .attr("height", y(20) - y(30))
      .attr("fill", "#ffffbf")
      .attr("opacity", 0.3);
    
    svg.append("rect")
      .attr("x", 0)
      .attr("y", y(20))
      .attr("width", width)
      .attr("height", y(0) - y(20))
      .attr("fill", "#fdae61")
      .attr("opacity", 0.3);
    
    // Add quality zone labels
    svg.append("text")
      .attr("x", 5)
      .attr("y", y(35))
      .attr("fill", "#1a9641")
      .attr("font-size", "10px")
      .text("Good quality");
    
    svg.append("text")
      .attr("x", 5)
      .attr("y", y(25))
      .attr("fill", "#a6611a")
      .attr("font-size", "10px")
      .text("Reasonable quality");
    
    svg.append("text")
      .attr("x", 5)
      .attr("y", y(10))
      .attr("fill", "#d7191c")
      .attr("font-size", "10px")
      .text("Poor quality");
    
    // Add boxplot for each position
    svg.selectAll("g.boxplot")
      .data(stats.positionData)
      .join("g")
        .attr("class", "boxplot")
        .attr("transform", d => `translate(${x(d.position)}, 0)`)
        .each(function(d) {
          const g = d3.select(this);
          
          // Box
          g.append("rect")
            .attr("x", -3)
            .attr("y", y(d.q3))
            .attr("width", 6)
            .attr("height", y(d.q1) - y(d.q3))
            .attr("fill", "#69b3a2")
            .attr("stroke", "#444444");
          
          // Median
          g.append("line")
            .attr("x1", -3)
            .attr("x2", 3)
            .attr("y1", y(d.median))
            .attr("y2", y(d.median))
            .attr("stroke", "#444444")
            .attr("stroke-width", 2);
          
          // Min line
          g.append("line")
            .attr("x1", 0)
            .attr("x2", 0)
            .attr("y1", y(d.q1))
            .attr("y2", y(d.min))
            .attr("stroke", "#444444")
            .attr("stroke-width", 1);
          
          // Max line
          g.append("line")
            .attr("x1", 0)
            .attr("x2", 0)
            .attr("y1", y(d.q3))
            .attr("y2", y(d.max))
            .attr("stroke", "#444444")
            .attr("stroke-width", 1);
          
          // Min cap
          g.append("line")
            .attr("x1", -3)
            .attr("x2", 3)
            .attr("y1", y(d.min))
            .attr("y2", y(d.min))
            .attr("stroke", "#444444")
            .attr("stroke-width", 1);
          
          // Max cap
          g.append("line")
            .attr("x1", -3)
            .attr("x2", 3)
            .attr("y1", y(d.max))
            .attr("y2", y(d.max))
            .attr("stroke", "#444444")
            .attr("stroke-width", 1);
        });
    
    // Add mean line
    svg.append("path")
      .datum(stats.positionData)
      .attr("fill", "none")
      .attr("stroke", "red")
      .attr("stroke-width", 2)
      .attr("d", d3.line()
        .x(d => x(d.position))
        .y(d => y(d.mean))
      );
    
    // Add title
    svg.append("text")
      .attr("x", width / 2)
      .attr("y", -10)
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .text("Quality Scores by Position");
    
  }, [stats]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Stats summary */}
      {stats && (
        <div className="border-b p-2 grid grid-cols-4 gap-4 text-sm">
          <div>
            <div className="font-medium">Read Count</div>
            <div>{stats.readCount}</div>
          </div>
          <div>
            <div className="font-medium">Max Length</div>
            <div>{stats.maxReadLength}</div>
          </div>
          <div>
            <div className="font-medium">Avg Quality</div>
            <div>{stats.overallStats.mean.toFixed(1)}</div>
          </div>
          <div>
            <div className="font-medium">Quality Range</div>
            <div>{stats.overallStats.min} - {stats.overallStats.max}</div>
          </div>
        </div>
      )}
      
      {/* Plot */}
      <div className="flex-1 p-2">
        {reads.length > 0 ? (
          <div className="h-full" ref={chartRef}></div>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No FASTQ data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Feature map viewer for GenBank, EMBL, GFF, BED files
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const FeatureMapViewer = ({ content, data, format, settings }) => {
  const chartRef = useRef(null);
  const [features, setFeatures] = useState([]);
  const [sequenceLength, setSequenceLength] = useState(0);
  const [zoomRange, setZoomRange] = useState({ start: 0, end: 0 });
  const [brushSelection, setBrushSelection] = useState(null);
  const [hoveredFeature, setHoveredFeature] = useState(null);
  
  useEffect(() => {
    // Extract features based on format
    if (content) {
      let extractedFeatures = [];
      let length = 0;
      
      if (format === 'genbank') {
        const result = extractGenBankFeatures(content);
        extractedFeatures = result.features;
        length = result.length;
      } else if (format === 'embl') {
        const result = extractEMBLFeatures(content);
        extractedFeatures = result.features;
        length = result.length;
      } else if (format === 'gff') {
        const result = extractGFFFeatures(content);
        extractedFeatures = result.features;
        length = result.length;
      } else if (format === 'bed') {
        const result = extractBEDFeatures(content);
        extractedFeatures = result.features;
        length = result.length;
      }
      
      setFeatures(extractedFeatures);
      setSequenceLength(length);
      setZoomRange({ start: 1, end: length });
    }
  }, [content, format]);
  
  // Extract features from GenBank
  const extractGenBankFeatures = (content) => {
    const features = [];
    let sequenceLength = 0;
    
    // Extract sequence length
    const locusMatch = content.match(/LOCUS\s+(\S+)\s+(\d+)\s+bp/);
    if (locusMatch) {
      sequenceLength = parseInt(locusMatch[2]);
    }
    
    // Extract features section
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
          
          // Parse location
          let start = 1;
          let end = sequenceLength;
          let strand = '+';
          
          // Simple range (e.g., "100..200")
          const simpleRangeMatch = location.match(/^(\d+)\.\.(\d+)$/);
          if (simpleRangeMatch) {
            start = parseInt(simpleRangeMatch[1]);
            end = parseInt(simpleRangeMatch[2]);
          }
          
          // Complement (e.g., "complement(100..200)")
          const complementMatch = location.match(/complement\((\d+)\.\.(\d+)\)/);
          if (complementMatch) {
            start = parseInt(complementMatch[1]);
            end = parseInt(complementMatch[2]);
            strand = '-';
          }
          
          // Single position (e.g., "42")
          const singlePosMatch = location.match(/^(\d+)$/);
          if (singlePosMatch) {
            start = parseInt(singlePosMatch[1]);
            end = start;
          }
          
          currentFeature = { 
            type, 
            start, 
            end, 
            strand,
            length: end - start + 1,
            qualifiers: {} 
          };
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
    
    return { features, length: sequenceLength };
  };
  
  // Extract features from EMBL
  const extractEMBLFeatures = (content) => {
    const features = [];
    let sequenceLength = 0;
    
    // Extract sequence length
    const idMatch = content.match(/ID.*?;\s+(\d+)\s+BP/);
    if (idMatch) {
      sequenceLength = parseInt(idMatch[1]);
    }
    
    // Extract features section
    const featureMatch = content.match(/FH.*?\nFT([\s\S]*?)(?=XX|$)/s);
    if (featureMatch) {
      const featureText = featureMatch[1];
      const featureLines = featureText.split('\n');
      
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
            
            // Parse location
            let start = 1;
            let end = sequenceLength;
            let strand = '+';
            
            // Simple range (e.g., "100..200")
            const simpleRangeMatch = location.match(/^(\d+)\.\.(\d+)$/);
            if (simpleRangeMatch) {
              start = parseInt(simpleRangeMatch[1]);
              end = parseInt(simpleRangeMatch[2]);
            }
            
            // Complement (e.g., "complement(100..200)")
            const complementMatch = location.match(/complement\((\d+)\.\.(\d+)\)/);
            if (complementMatch) {
              start = parseInt(complementMatch[1]);
              end = parseInt(complementMatch[2]);
              strand = '-';
            }
            
            currentFeature = { 
              type, 
              start, 
              end, 
              strand,
              length: end - start + 1,
              qualifiers: {} 
            };
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
    
    return { features, length: sequenceLength };
  };
  
  // Extract features from GFF
  const extractGFFFeatures = (content) => {
    const features = [];
    let sequenceLength = 0;
    
    // Look for sequence-region pragma to get sequence length
    const regionMatch = content.match(/##sequence-region\s+\S+\s+\d+\s+(\d+)/);
    if (regionMatch) {
      sequenceLength = parseInt(regionMatch[1]);
    }
    
    // Parse GFF lines
    const lines = content.split('\n');
    
    for (const line of lines) {
      // Skip comments and empty lines
      if (line.startsWith('#') || !line.trim()) {
        continue;
      }
      
      const fields = line.split('\t');
      if (fields.length < 8) {
        continue;
      }
      
      const seqId = fields[0];
      const source = fields[1];
      const type = fields[2];
      const start = parseInt(fields[3]);
      const end = parseInt(fields[4]);
      const score = fields[5];
      const strand = fields[6];
      const phase = fields[7];
      const attributes = fields[8] || '';
      
      // Update sequence length if needed
      if (end > sequenceLength) {
        sequenceLength = end;
      }
      
      // Parse attributes
      const attrMap = {};
      const attrPairs = attributes.split(';');
      
      for (const pair of attrPairs) {
        const [key, value] = pair.split('=');
        if (key && value) {
          attrMap[key] = value;
        }
      }
      
      // Create feature
      features.push({
        type,
        start,
        end,
        strand,
        length: end - start + 1,
        source,
        score: score !== '.' ? parseFloat(score) : null,
        qualifiers: attrMap
      });
    }
    
    return { features, length: sequenceLength };
  };
  
  // Extract features from BED
  const extractBEDFeatures = (content) => {
    const features = [];
    let sequenceLength = 0;
    
    // Parse BED lines
    const lines = content.split('\n');
    
    for (const line of lines) {
      // Skip comments and empty lines
      if (line.startsWith('#') || !line.trim()) {
        continue;
      }
      
      const fields = line.split('\t');
      if (fields.length < 3) {
        continue;
      }
      
      const chrom = fields[0];
      const start = parseInt(fields[1]); // BED is 0-based
      const end = parseInt(fields[2]);
      const name = fields.length > 3 ? fields[3] : '';
      const score = fields.length > 4 ? fields[4] : '.';
      const strand = fields.length > 5 ? fields[5] : '+';
      
      // Update sequence length if needed
      if (end > sequenceLength) {
        sequenceLength = end;
      }
      
      // Create feature
      features.push({
        type: name || 'feature',
        start: start + 1, // Convert to 1-based
        end,
        strand,
        length: end - start,
        chrom,
        score: score !== '.' ? parseFloat(score) : null,
        qualifiers: {}
      });
    }
    
    return { features, length: sequenceLength };
  };
  
  // Get color for feature type
  const getFeatureColor = (type) => {
    const colorMap = {
      'gene': '#1f77b4',
      'CDS': '#ff7f0e',
      'mRNA': '#2ca02c',
      'exon': '#d62728',
      'tRNA': '#9467bd',
      'rRNA': '#8c564b',
      'misc_feature': '#e377c2',
      'regulatory': '#7f7f7f',
      'source': '#bcbd22',
      'repeat_region': '#17becf'
    };
    
    // Return color or default
    return colorMap[type] || '#aaaaaa';
  };
  
  // Render feature map
  useEffect(() => {
    if (!chartRef.current || features.length === 0 || sequenceLength === 0) return;
    
    // Clear previous chart
    d3.select(chartRef.current).selectAll("*").remove();
    
    // Set up dimensions
    const margin = { top: 40, right: 30, bottom: 60, left: 100 };
    const width = chartRef.current.clientWidth - margin.left - margin.right;
    const height = chartRef.current.clientHeight - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(chartRef.current)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // X scale
    const x = d3.scaleLinear()
      .domain([zoomRange.start, zoomRange.end])
      .range([0, width]);
    
    // Group features by type
    const featureTypes = [...new Set(features.map(f => f.type))];
    
    // Y scale
    const y = d3.scaleBand()
      .domain(featureTypes)
      .range([0, height])
      .padding(0.3);
    
    // X axis
    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(x).ticks(10))
      .append("text")
        .attr("x", width / 2)
        .attr("y", margin.bottom - 10)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Position");
    
    // Y axis
    svg.append("g")
      .call(d3.axisLeft(y))
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left + 15)
        .attr("x", -height / 2)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Feature Type");
    
    // Tooltip
    const tooltip = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0)
      .style("background-color", "white")
      .style("border", "solid")
      .style("border-width", "1px")
      .style("border-radius", "5px")
      .style("padding", "10px")
      .style("position", "absolute");
    
    // Add features
    svg.selectAll("rect.feature")
      .data(features)
      .join("rect")
        .attr("class", "feature")
        .attr("x", d => Math.max(0, x(d.start)))
        .attr("y", d => y(d.type))
        .attr("width", d => Math.max(1, x(d.end) - x(d.start)))
        .attr("height", y.bandwidth())
        .attr("fill", d => getFeatureColor(d.type))
        .attr("stroke", "black")
        .attr("stroke-width", 0.5)
        .attr("rx", 3) // Rounded corners
        .on("mouseover", function(d) {
          // Highlight feature
          d3.select(this)
            .attr("stroke-width", 2)
            .attr("stroke", "#ff0000");
          
          // Show tooltip
          tooltip.transition()
            .duration(200)
            .style("opacity", .9);
          
          let tooltipContent = `<strong>${d.type}</strong><br>`;
          tooltipContent += `Position: ${d.start} - ${d.end}<br>`;
          tooltipContent += `Length: ${d.length} bp<br>`;
          tooltipContent += `Strand: ${d.strand}<br>`;
          
          // Add qualifiers
          if (Object.keys(d.qualifiers).length > 0) {
            tooltipContent += "<hr>";
            for (const [key, value] of Object.entries(d.qualifiers)) {
              tooltipContent += `<strong>${key}:</strong> ${value}<br>`;
            }
          }
          
          tooltip.html(tooltipContent)
            .style("left", (event.pageX + 10) + "px")
            .style("top", (event.pageY - 28) + "px");
          
          setHoveredFeature(d);
        })
        .on("mouseout", function(d) {
          // Remove highlight
          d3.select(this)
            .attr("stroke-width", 0.5)
            .attr("stroke", "black");
          
          // Hide tooltip
          tooltip.transition()
            .duration(500)
            .style("opacity", 0);
          
          setHoveredFeature(null);
        });
    
    // Add strand markers
    svg.selectAll("path.strand")
      .data(features.filter(f => f.strand === '-'))
      .join("path")
        .attr("class", "strand")
        .attr("d", d => {
          const xPos = Math.max(0, x(d.start));
          const yPos = y(d.type) + y.bandwidth() / 2;
          const markerSize = Math.min(y.bandwidth() / 2, 8);
          return `M ${xPos + markerSize * 2} ${yPos} L ${xPos} ${yPos - markerSize} L ${xPos} ${yPos + markerSize} Z`;
        })
        .attr("fill", "black");
    
    // Add title
    svg.append("text")
      .attr("x", width / 2)
      .attr("y", -20)
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .text("Feature Map");
    
    // Zoom brush
    const brush = d3.brushX()
      .extent([[0, -margin.top / 2], [width, 0]])
      .on("brush", function(event) {
        setBrushSelection(event.selection);
      })
      .on("end", function(event) {
        if (!event.selection) return;
        const [x0, x1] = event.selection.map(x.invert);
        setZoomRange({ start: Math.max(1, Math.floor(x0)), end: Math.ceil(x1) });
        setBrushSelection(null);
      });
    
    // Add context view (mini-map)
    const contextHeight = 20;
    const contextY = -margin.top / 2 - contextHeight;
    
    // Context X scale
    const xContext = d3.scaleLinear()
      .domain([1, sequenceLength])
      .range([0, width]);
    
    // Add context background
    svg.append("rect")
      .attr("x", 0)
      .attr("y", contextY)
      .attr("width", width)
      .attr("height", contextHeight)
      .attr("fill", "#f0f0f0");
    
    // Add context features
    svg.selectAll("rect.context-feature")
      .data(features)
      .join("rect")
        .attr("class", "context-feature")
        .attr("x", d => xContext(d.start))
        .attr("y", contextY)
        .attr("width", d => Math.max(1, xContext(d.end) - xContext(d.start)))
        .attr("height", contextHeight)
        .attr("fill", d => getFeatureColor(d.type))
        .attr("opacity", 0.7);
    
    // Add context axis
    svg.append("g")
      .attr("transform", `translate(0,${contextY + contextHeight})`)
      .call(d3.axisBottom(xContext).ticks(5));
    
    // Add brush
    svg.append("g")
      .attr("class", "brush")
      .attr("transform", `translate(0,${contextY})`)
      .call(brush)
      .call(brush.move, [xContext(zoomRange.start), xContext(zoomRange.end)]);
    
    // Add brush handle styles
    svg.selectAll(".brush>.handle")
      .attr("fill", "#69b3a2")
      .attr("stroke", "#000000");
    
    // Return cleanup function
    return () => {
      d3.select("body").selectAll(".tooltip").remove();
    };
  }, [features, sequenceLength, zoomRange]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div>
          <span className="text-sm text-gray-600 mr-2">Features: {features.length}</span>
          <span className="text-sm text-gray-600">Sequence Length: {sequenceLength} bp</span>
        </div>
        
        <div className="flex items-center space-x-2">
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm"
            onClick={() => setZoomRange({ start: 1, end: sequenceLength })}
          >
            Reset Zoom
          </button>
        </div>
      </div>
      
      {/* Feature details panel for hovered feature */}
      {hoveredFeature && (
        <div className="border-b p-2 bg-blue-50">
          <h3 className="font-medium">{hoveredFeature.type}</h3>
          <div className="grid grid-cols-3 text-sm">
            <div>Position: {hoveredFeature.start} - {hoveredFeature.end}</div>
            <div>Length: {hoveredFeature.length} bp</div>
            <div>Strand: {hoveredFeature.strand}</div>
          </div>
          {Object.keys(hoveredFeature.qualifiers).length > 0 && (
            <div className="mt-1 text-sm">
              <div className="font-medium">Qualifiers:</div>
              <div className="grid grid-cols-2 gap-1">
                {Object.entries(hoveredFeature.qualifiers).map(([key, value]) => (
                  <div key={key} className="overflow-hidden text-ellipsis whitespace-nowrap">
                    <span className="font-medium">{key}:</span> {value}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
      
      {/* Feature map */}
      <div className="flex-1 p-2">
        {features.length > 0 ? (
          <div className="h-full" ref={chartRef}></div>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No feature data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Circular genome viewer component
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const CircularViewer = ({ content, data, settings }) => {
  const svgRef = useRef(null);
  const [features, setFeatures] = useState([]);
  const [sequenceLength, setSequenceLength] = useState(0);
  const [trackCount, setTrackCount] = useState(1);
  const [rotation, setRotation] = useState(0);
  const [hoveredFeature, setHoveredFeature] = useState(null);
  const [selectedFeatureTypes, setSelectedFeatureTypes] = useState([]);
  
  useEffect(() => {
    // Extract features based on format
    if (content) {
      let extractedFeatures = [];
      let length = 0;
      
      if (data?.format === 'genbank' || data?.format === 'embl') {
        // Use the features from the data if available
        if (data.features && data.sequenceLength) {
          extractedFeatures = data.features;
          length = data.sequenceLength;
        } else {
          // Otherwise extract them from content
          const result = extractAnnotatedFeatures(content, data.format);
          extractedFeatures = result.features;
          length = result.length;
        }
      }
      
      setFeatures(extractedFeatures);
      setSequenceLength(length);
      
      // Get unique feature types for filtering
      const types = [...new Set(extractedFeatures.map(f => f.type))];
      setSelectedFeatureTypes(types);
    }
  }, [content, data]);
  
  // Extract features from GenBank/EMBL
  const extractAnnotatedFeatures = (content, format) => {
    const features = [];
    let sequenceLength = 0;
    
    if (format === 'genbank') {
      // Extract sequence length
      const locusMatch = content.match(/LOCUS\s+(\S+)\s+(\d+)\s+bp/);
      if (locusMatch) {
        sequenceLength = parseInt(locusMatch[2]);
      }
      
      // Extract features section
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
            
            // Parse location
            let start = 1;
            let end = sequenceLength;
            let strand = '+';
            
            // Simple range (e.g., "100..200")
            const simpleRangeMatch = location.match(/^(\d+)\.\.(\d+)$/);
            if (simpleRangeMatch) {
              start = parseInt(simpleRangeMatch[1]);
              end = parseInt(simpleRangeMatch[2]);
            }
            
            // Complement (e.g., "complement(100..200)")
            const complementMatch = location.match(/complement\((\d+)\.\.(\d+)\)/);
            if (complementMatch) {
              start = parseInt(complementMatch[1]);
              end = parseInt(complementMatch[2]);
              strand = '-';
            }
            
            currentFeature = { 
              type, 
              start, 
              end, 
              strand,
              length: end - start + 1,
              qualifiers: {} 
            };
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
      // Similar parsing for EMBL format
      // ...
    }
    
    return { features, length: sequenceLength };
  };
  
  // Get color for feature type
  const getFeatureColor = (type) => {
    const colorMap = {
      'gene': '#1f77b4',
      'CDS': '#ff7f0e',
      'mRNA': '#2ca02c',
      'exon': '#d62728',
      'tRNA': '#9467bd',
      'rRNA': '#8c564b',
      'misc_feature': '#e377c2',
      'regulatory': '#7f7f7f',
      'source': '#bcbd22',
      'repeat_region': '#17becf'
    };
    
    // Return color or default
    return colorMap[type] || '#aaaaaa';
  };
  
  // Render circular map
  useEffect(() => {
    if (!svgRef.current || features.length === 0 || sequenceLength === 0) return;
    
    // Clear previous content
    d3.select(svgRef.current).selectAll("*").remove();
    
    // Filter features by selected types
    const filteredFeatures = features.filter(f => selectedFeatureTypes.includes(f.type));
    
    // Set up dimensions
    const svg = d3.select(svgRef.current);
    const width = svgRef.current.clientWidth;
    const height = svgRef.current.clientHeight;
    const margin = 60;
    const radius = Math.min(width, height) / 2 - margin;
    
    // Create a group at the center
    const g = svg.append("g")
      .attr("transform", `translate(${width / 2},${height / 2}) rotate(${rotation})`);
    
    // Add genome circle
    g.append("circle")
      .attr("r", radius)
      .attr("fill", "none")
      .attr("stroke", "#cccccc")
      .attr("stroke-width", 2);
    
    // Radial scale
    const angle = d3.scaleLinear()
      .domain([0, sequenceLength])
      .range([0, 2 * Math.PI]);
    
    // Group features by type to determine tracks
    const featuresByType = d3.group(filteredFeatures, d => d.type);
    const types = Array.from(featuresByType.keys());
    
    // Determine radial positions for tracks
    const trackWidth = 20;
    const maxTracks = Math.floor((radius - 20) / trackWidth);
    const trackRadius = {};
    
    // Assign tracks by feature type
    types.forEach((type, i) => {
      trackRadius[type] = radius - (i % maxTracks + 1) * trackWidth;
    });
    
    // Track count
    setTrackCount(Math.min(types.length, maxTracks));
    
    // Add position ticks
    const tickStep = sequenceLength > 10000 ? 1000 : 100;
    const ticks = d3.range(0, sequenceLength, tickStep);
    
    g.selectAll("line.tick")
      .data(ticks)
      .enter()
      .append("line")
        .attr("class", "tick")
        .attr("x1", 0)
        .attr("y1", 0)
        .attr("x2", (d) => (radius + 10) * Math.cos(angle(d) - Math.PI / 2))
        .attr("y2", (d) => (radius + 10) * Math.sin(angle(d) - Math.PI / 2))
        .attr("stroke", "#999999")
        .attr("stroke-width", d => d % (tickStep * 10) === 0 ? 2 : 1);
    
    // Add position labels
    g.selectAll("text.tick-label")
      .data(ticks.filter(d => d % (tickStep * 10) === 0))
      .enter()
      .append("text")
        .attr("class", "tick-label")
        .attr("x", (d) => (radius + 25) * Math.cos(angle(d) - Math.PI / 2))
        .attr("y", (d) => (radius + 25) * Math.sin(angle(d) - Math.PI / 2))
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("font-size", "10px")
        .text(d => d.toLocaleString());
    
    // Add features
    types.forEach(type => {
      const typeFeatures = featuresByType.get(type);
      const trackWidth = 10;
      const currentRadius = trackRadius[type];
      
      // Draw arcs for each feature
      g.selectAll(`path.feature-${type}`)
        .data(typeFeatures)
        .enter()
        .append("path")
        .attr("class", `feature-${type}`)
        .attr("d", d => {
          const startAngle = angle(d.start);
          const endAngle = angle(d.end);
          
          const arc = d3.arc()
            .innerRadius(currentRadius - trackWidth / 2)
            .outerRadius(currentRadius + trackWidth / 2)
            .startAngle(startAngle - Math.PI / 2)
            .endAngle(endAngle - Math.PI / 2);
          
          return arc();
        })
        .attr("fill", getFeatureColor(type))
        .attr("stroke", "white")
        .attr("stroke-width", 0.5)
        .on("mouseover", function(d) {
          // Highlight
          d3.select(this)
            .attr("stroke", "black")
            .attr("stroke-width", 2);
          
          setHoveredFeature(d);
        })
        .on("mouseout", function() {
          // Remove highlight
          d3.select(this)
            .attr("stroke", "white")
            .attr("stroke-width", 0.5);
          
          setHoveredFeature(null);
        });
      
      // Add feature type label
      const labelAngle = -Math.PI / 2; // Top position
      g.append("text")
        .attr("x", (currentRadius) * Math.cos(labelAngle))
        .attr("y", (currentRadius) * Math.sin(labelAngle))
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("font-size", "8px")
        .attr("font-weight", "bold")
        .attr("transform", `rotate(90, ${(currentRadius) * Math.cos(labelAngle)}, ${(currentRadius) * Math.sin(labelAngle)})`)
        .text(type);
    });
    
    // Add genome name and length info
    g.append("text")
      .attr("x", 0)
      .attr("y", 0)
      .attr("text-anchor", "middle")
      .attr("alignment-baseline", "middle")
      .attr("font-weight", "bold")
      .text(`${sequenceLength.toLocaleString()} bp`);
    
  }, [features, sequenceLength, selectedFeatureTypes, rotation]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div>
          <span className="text-sm text-gray-600 mr-2">
            Features: {features.filter(f => selectedFeatureTypes.includes(f.type)).length}
          </span>
          <span className="text-sm text-gray-600">Tracks: {trackCount}</span>
        </div>
        
        <div className="flex items-center space-x-2">
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm"
            onClick={() => setRotation((rotation + 30) % 360)}
          >
            Rotate
          </button>
          
          <select
            className="text-sm border rounded p-1"
            value={selectedFeatureTypes.length === 0 ? "" : "custom"}
            onChange={(e) => {
              if (e.target.value === "all") {
                setSelectedFeatureTypes([...new Set(features.map(f => f.type))]);
              } else if (e.target.value === "none") {
                setSelectedFeatureTypes([]);
              }
            }}
          >
            <option value="custom">
              {selectedFeatureTypes.length} types selected
            </option>
            <option value="all">Select all types</option>
            <option value="none">Clear selection</option>
          </select>
        </div>
      </div>
      
      {/* Feature type selection */}
      <div className="border-b p-2 flex flex-wrap gap-1">
        {[...new Set(features.map(f => f.type))].map(type => (
          <label key={type} className="inline-flex items-center cursor-pointer">
            <input
              type="checkbox"
              className="form-checkbox h-3 w-3"
              checked={selectedFeatureTypes.includes(type)}
              onChange={(e) => {
                if (e.target.checked) {
                  setSelectedFeatureTypes([...selectedFeatureTypes, type]);
                } else {
                  setSelectedFeatureTypes(selectedFeatureTypes.filter(t => t !== type));
                }
              }}
            />
            <span className="ml-1 mr-2 text-xs">
              <span className="inline-block w-2 h-2 mr-1" style={{backgroundColor: getFeatureColor(type)}}></span>
              {type}
            </span>
          </label>
        ))}
      </div>
      
      {/* Feature details panel for hovered feature */}
      {hoveredFeature && (
        <div className="border-b p-2 bg-blue-50">
          <h3 className="font-medium">{hoveredFeature.type}</h3>
          <div className="grid grid-cols-3 text-sm">
            <div>Position: {hoveredFeature.start} - {hoveredFeature.end}</div>
            <div>Length: {hoveredFeature.length} bp</div>
            <div>Strand: {hoveredFeature.strand}</div>
          </div>
          {Object.keys(hoveredFeature.qualifiers).length > 0 && (
            <div className="mt-1 text-sm">
              <div className="font-medium">Qualifiers:</div>
              <div className="grid grid-cols-2 gap-1">
                {Object.entries(hoveredFeature.qualifiers).map(([key, value]) => (
                  <div key={key} className="overflow-hidden text-ellipsis whitespace-nowrap">
                    <span className="font-medium">{key}:</span> {value}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
      
      {/* SVG container */}
      <div className="flex-1 overflow-hidden">
        {features.length > 0 ? (
          <svg ref={svgRef} width="100%" height="100%">
            {/* SVG content will be added by D3 */}
          </svg>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No feature data available for circular view</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * 3D Structure viewer for PDB/mmCIF files
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const StructureViewer = ({ content, data, settings }) => {
  const containerRef = useRef(null);
  const [loading, setLoading] = useState(true);
  const [chains, setChains] = useState([]);
  const [displayMode, setDisplayMode] = useState('cartoon');
  const [colorMode, setColorMode] = useState('chain');
  
  useEffect(() => {
    if (!containerRef.current || !content) return;
    
    // Initialize 3D renderer
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0xffffff);
    
    const camera = new THREE.PerspectiveCamera(
      75,
      containerRef.current.clientWidth / containerRef.current.clientHeight,
      0.1,
      1000
    );
    camera.position.z = 50;
    
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(containerRef.current.clientWidth, containerRef.current.clientHeight);
    containerRef.current.appendChild(renderer.domElement);
    
    // Add light
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);
    
    const light = new THREE.DirectionalLight(0xffffff, 1);
    light.position.set(10, 10, 10);
    scene.add(light);
    
    // Add protein structure (simplified)
    // In a real implementation, we'd need a full PDB parser
    // Here's a placeholder visualization with a sphere model
    
    // Extract atom coordinates from PDB
    const atoms = [];
    const chainsSet = new Set();
    
    if (data?.format === 'pdb') {
      const lines = content.split('\n');
      
      for (const line of lines) {
        if (line.startsWith('ATOM')) {
          const x = parseFloat(line.substring(30, 38));
          const y = parseFloat(line.substring(38, 46));
          const z = parseFloat(line.substring(46, 54));
          const atomName = line.substring(12, 16).trim();
          const residue = line.substring(17, 20).trim();
          const chainId = line.substring(21, 22).trim();
          const residueNum = parseInt(line.substring(22, 26));
          
          // Add atom
          atoms.push({
            x, y, z, atomName, residue, chainId, residueNum
          });
          
          // Track chains
          chainsSet.add(chainId);
        }
      }
    }
    
    // Update chains state
    setChains(Array.from(chainsSet));
    
    // Compute center of geometry
    let center = { x: 0, y: 0, z: 0 };
    if (atoms.length > 0) {
      for (const atom of atoms) {
        center.x += atom.x;
        center.y += atom.y;
        center.z += atom.z;
      }
      center.x /= atoms.length;
      center.y /= atoms.length;
      center.z /= atoms.length;
    }
    
    // Function to get color based on chain
    const getChainColor = (chainId) => {
      const colors = [
        0x1f77b4, 0xff7f0e, 0x2ca02c, 0xd62728,
        0x9467bd, 0x8c564b, 0xe377c2, 0x7f7f7f
      ];
      
      // Hash the chain ID to a consistent color
      const index = chainId.charCodeAt(0) % colors.length;
      return colors[index];
    };
    
    // Add atoms as spheres
    const geometry = new THREE.SphereGeometry(0.4, 8, 8);
    
    for (const atom of atoms) {
      const color = getChainColor(atom.chainId);
      const material = new THREE.MeshLambertMaterial({ color });
      
      const sphere = new THREE.Mesh(geometry, material);
      sphere.position.set(
        atom.x - center.x,
        atom.y - center.y,
        atom.z - center.z
      );
      
      scene.add(sphere);
    }
    
    // Add simple animation
    const animate = () => {
      requestAnimationFrame(animate);
      
      scene.rotation.y += 0.005;
      
      renderer.render(scene, camera);
    };
    
    animate();
    setLoading(false);
    
    // Cleanup
    return () => {
      if (containerRef.current && containerRef.current.contains(renderer.domElement)) {
        containerRef.current.removeChild(renderer.domElement);
      }
      renderer.dispose();
    };
  }, [content, data, displayMode, colorMode]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 grid grid-cols-3 gap-2">
        <div>
          <label className="block text-xs text-gray-600 mb-1">Display Style</label>
          <select
            className="w-full text-sm p-1 border rounded"
            value={displayMode}
            onChange={(e) => setDisplayMode(e.target.value)}
          >
            <option value="cartoon">Cartoon</option>
            <option value="sphere">Sphere</option>
            <option value="stick">Stick</option>
            <option value="line">Line</option>
            <option value="ribbon">Ribbon</option>
          </select>
        </div>
        
        <div>
          <label className="block text-xs text-gray-600 mb-1">Color By</label>
          <select
            className="w-full text-sm p-1 border rounded"
            value={colorMode}
            onChange={(e) => setColorMode(e.target.value)}
          >
            <option value="chain">Chain</option>
            <option value="residue">Residue Type</option>
            <option value="element">Element</option>
            <option value="bfactor">B-Factor</option>
            <option value="hydrophobicity">Hydrophobicity</option>
          </select>
        </div>
        
        <div>
          <label className="block text-xs text-gray-600 mb-1">Chains</label>
          <select
            className="w-full text-sm p-1 border rounded"
            defaultValue="all"
          >
            <option value="all">All Chains</option>
            {chains.map(chain => (
              <option key={chain} value={chain}>Chain {chain}</option>
            ))}
          </select>
        </div>
      </div>
      
      {/* 3D Viewer */}
      <div className="flex-1" ref={containerRef}>
        {loading && (
          <div className="flex items-center justify-center h-full">
            <RefreshCw className="w-6 h-6 text-blue-600 animate-spin" />
            <span className="ml-2 text-gray-600">Loading 3D structure...</span>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Contact map viewer for protein structures
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const ContactMapViewer = ({ content, data, settings }) => {
  const svgRef = useRef(null);
  const [atoms, setAtoms] = useState([]);
  const [contacts, setContacts] = useState([]);
  const [threshold, setThreshold] = useState(8.0); // Angstroms
  const [selectedChain, setSelectedChain] = useState('all');
  const [chains, setChains] = useState([]);
  
  useEffect(() => {
    if (!content) return;
    
    // Extract atoms from PDB
    const extractedAtoms = [];
    const chainsSet = new Set();
    
    if (data?.format === 'pdb') {
      const lines = content.split('\n');
      
      for (const line of lines) {
        if (line.startsWith('ATOM')) {
          const x = parseFloat(line.substring(30, 38));
          const y = parseFloat(line.substring(38, 46));
          const z = parseFloat(line.substring(46, 54));
          const atomName = line.substring(12, 16).trim();
          const residue = line.substring(17, 20).trim();
          const chainId = line.substring(21, 22).trim();
          const residueNum = parseInt(line.substring(22, 26));
          
          // Only add CA atoms for simplicity
          if (atomName === 'CA') {
            extractedAtoms.push({
              x, y, z, atomName, residue, chainId, residueNum
            });
            
            // Track chains
            chainsSet.add(chainId);
          }
        }
      }
    }
    
    setAtoms(extractedAtoms);
    setChains(Array.from(chainsSet));
    
    // Calculate contacts
    calculateContacts(extractedAtoms, threshold);
  }, [content, data]);
  
  // Re-calculate contacts when threshold changes
  useEffect(() => {
    calculateContacts(atoms, threshold);
  }, [atoms, threshold]);
  
  // Calculate residue-residue contacts
  const calculateContacts = (atoms, threshold) => {
    const contacts = [];
    
    // Filter by selected chain if needed
    const filteredAtoms = selectedChain === 'all' 
      ? atoms 
      : atoms.filter(atom => atom.chainId === selectedChain);
    
    // Calculate distances between residues
    for (let i = 0; i < filteredAtoms.length; i++) {
      for (let j = i + 1; j < filteredAtoms.length; j++) {
        const atom1 = filteredAtoms[i];
        const atom2 = filteredAtoms[j];
        
        // Skip nearby residues in sequence (trivial contacts)
        if (atom1.chainId === atom2.chainId && 
            Math.abs(atom1.residueNum - atom2.residueNum) < 4) {
          continue;
        }
        
        // Calculate Euclidean distance
        const distance = Math.sqrt(
          Math.pow(atom1.x - atom2.x, 2) +
          Math.pow(atom1.y - atom2.y, 2) +
          Math.pow(atom1.z - atom2.z, 2)
        );
        
        // Check if below threshold
        if (distance <= threshold) {
          contacts.push({
            residue1: atom1.residueNum,
            chain1: atom1.chainId,
            residue2: atom2.residueNum,
            chain2: atom2.chainId,
            distance
          });
        }
      }
    }
    
    setContacts(contacts);
  };
  
  // Render contact map
  useEffect(() => {
    if (!svgRef.current || atoms.length === 0) return;
    
    // Clear previous map
    d3.select(svgRef.current).selectAll("*").remove();
    
    // Set up dimensions
    const margin = { top: 50, right: 50, bottom: 70, left: 70 };
    const width = svgRef.current.clientWidth - margin.left - margin.right;
    const height = svgRef.current.clientHeight - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(svgRef.current)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // Find min and max residue numbers
    const filteredAtoms = selectedChain === 'all' 
      ? atoms 
      : atoms.filter(atom => atom.chainId === selectedChain);
    
    if (filteredAtoms.length === 0) return;
    
    const minResidue = d3.min(filteredAtoms, d => d.residueNum);
    const maxResidue = d3.max(filteredAtoms, d => d.residueNum);
    
    // Scales
    const x = d3.scaleLinear()
      .domain([minResidue, maxResidue])
      .range([0, width]);
    
    const y = d3.scaleLinear()
      .domain([minResidue, maxResidue])
      .range([height, 0]);
    
    // Add axes
    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(x).ticks(10))
      .append("text")
        .attr("x", width / 2)
        .attr("y", 40)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Residue Number");
    
    svg.append("g")
      .call(d3.axisLeft(y).ticks(10))
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -40)
        .attr("x", -height / 2)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Residue Number");
    
    // Add contact points
    svg.selectAll("circle")
      .data(contacts)
      .enter()
      .append("circle")
        .attr("cx", d => x(d.residue1))
        .attr("cy", d => y(d.residue2))
        .attr("r", 2)
        .attr("fill", d => {
          // Color based on distance
          const intensity = 1 - (d.distance / threshold);
          return d3.interpolateBlues(intensity);
        })
        .attr("opacity", 0.7)
        .on("mouseover", function(d) {
          // Highlight
          d3.select(this)
            .attr("r", 4)
            .attr("stroke", "black");
          
          // Tooltip
          svg.append("text")
            .attr("id", "tooltip")
            .attr("x", x(d.residue1) + 5)
            .attr("y", y(d.residue2) - 5)
            .attr("fill", "black")
            .style("font-size", "10px")
            .text(`${d.chain1}:${d.residue1} - ${d.chain2}:${d.residue2} (${d.distance.toFixed(1)}Å)`);
        })
        .on("mouseout", function() {
          // Remove highlight
          d3.select(this)
            .attr("r", 2)
            .attr("stroke", "none");
          
          // Remove tooltip
          svg.select("#tooltip").remove();
        });
    
    // Symmetry - add the other half of the map
    svg.selectAll("circle.sym")
      .data(contacts)
      .enter()
      .append("circle")
        .attr("cx", d => x(d.residue2))
        .attr("cy", d => y(d.residue1))
        .attr("r", 2)
        .attr("fill", d => {
          // Color based on distance
          const intensity = 1 - (d.distance / threshold);
          return d3.interpolateBlues(intensity);
        })
        .attr("opacity", 0.7)
        .on("mouseover", function(d) {
          // Highlight
          d3.select(this)
            .attr("r", 4)
            .attr("stroke", "black");
          
          // Tooltip
          svg.append("text")
            .attr("id", "tooltip")
            .attr("x", x(d.residue2) + 5)
            .attr("y", y(d.residue1) - 5)
            .attr("fill", "black")
            .style("font-size", "10px")
            .text(`${d.chain2}:${d.residue2} - ${d.chain1}:${d.residue1} (${d.distance.toFixed(1)}Å)`);
        })
        .on("mouseout", function() {
          // Remove highlight
          d3.select(this)
            .attr("r", 2)
            .attr("stroke", "none");
          
          // Remove tooltip
          svg.select("#tooltip").remove();
        });
    
    // Add title
    svg.append("text")
      .attr("x", width / 2)
      .attr("y", -20)
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .text(`Contact Map (${threshold.toFixed(1)}Å threshold)`);
    
  }, [atoms, contacts, threshold, selectedChain]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex items-center justify-between">
        <div className="flex items-center">
          <label className="text-sm text-gray-600 mr-2">Distance Threshold:</label>
          <input
            type="range"
            min="4"
            max="20"
            step="0.5"
            value={threshold}
            onChange={(e) => setThreshold(parseFloat(e.target.value))}
            className="w-32"
          />
          <span className="ml-2 text-sm">{threshold.toFixed(1)} Å</span>
        </div>
        
        <div className="flex items-center">
          <label className="text-sm text-gray-600 mr-2">Chain:</label>
          <select
            className="text-sm border rounded p-1"
            value={selectedChain}
            onChange={(e) => setSelectedChain(e.target.value)}
          >
            <option value="all">All Chains</option>
            {chains.map(chain => (
              <option key={chain} value={chain}>Chain {chain}</option>
            ))}
          </select>
        </div>
      </div>
      
      {/* Stats */}
      <div className="border-b p-2 text-sm text-gray-600">
        {atoms.length > 0 ? (
          <>
            <span>Residues: {atoms.length}</span>
            <span className="mx-2">|</span>
            <span>Contacts: {contacts.length}</span>
            <span className="mx-2">|</span>
            <span>Contact Density: {(contacts.length / (atoms.length * (atoms.length - 1) / 2) * 100).toFixed(2)}%</span>
          </>
        ) : (
          <span>No structure data available</span>
        )}
      </div>
      
      {/* Map */}
      <div className="flex-1" ref={svgRef}>
        {atoms.length === 0 && (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No contact data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Alignment viewer for multiple sequence alignments
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const AlignmentViewer = ({ content, data, settings }) => {
  const [sequences, setSequences] = useState([]);
  const [conservation, setConservation] = useState('');
  const [zoomLevel, setZoomLevel] = useState(1);
  const [colorScheme, setColorScheme] = useState(settings.colorScheme || 'clustal');
  const [viewportStart, setViewportStart] = useState(0);
  const [viewportEnd, setViewportEnd] = useState(100);
  const [showingFullWidth, setShowingFullWidth] = useState(false);
  
  useEffect(() => {
    // Extract sequences based on format
    if (content) {
      let extractedSequences = [];
      let conservationStr = '';
      
      if (data?.format === 'clustal') {
        const result = extractClustalData(content);
        extractedSequences = result.sequences;
        conservationStr = result.conservation;
      } else if (data?.format === 'stockholm') {
        const result = extractStockholmData(content);
        extractedSequences = result.sequences;
        conservationStr = result.conservation;
      }
      
      setSequences(extractedSequences);
      setConservation(conservationStr);
      
      // Set viewport to show full width if alignment is small
      if (extractedSequences.length > 0) {
        const firstLen = extractedSequences[0].sequence.length;
        if (firstLen <= 100) {
          setViewportEnd(firstLen);
          setShowingFullWidth(true);
        } else {
          setViewportEnd(100);
          setShowingFullWidth(false);
        }
      }
    }
  }, [content, data]);
  
  // Extract data from CLUSTAL format
  const extractClustalData = (content) => {
    const sequences = [];
    let conservation = '';
    
    const lines = content.split('\n');
    
    // Find the header line
    const headerLineIndex = lines.findIndex(line => line.includes('CLUSTAL'));
    
    if (headerLineIndex !== -1) {
      const seqMap = {};
      const consLines = [];
      
      // Parse sequences and conservation lines
      for (let i = headerLineIndex + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        
        if (!line) continue;
        
        if (line.startsWith(' ')) {
          // Conservation line
          const consSymbols = line.trim();
          consLines.push(consSymbols);
        } else {
          // Sequence line
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
      }
      
      // Convert to array
      for (const [name, sequence] of Object.entries(seqMap)) {
        sequences.push({
          name,
          sequence,
          length: sequence.length
        });
      }
      
      // Join conservation lines
      conservation = consLines.join('');
    }
    
    return { sequences, conservation };
  };
  
  // Extract data from Stockholm format
  const extractStockholmData = (content) => {
    const sequences = [];
    let conservation = '';
    
    const lines = content.split('\n');
    
    const seqMap = {};
    
    for (const line of lines) {
      if (line.startsWith('#=GC RF')) {
        // Conservation line
        const parts = line.split(/\s+/);
        if (parts.length >= 3) {
          conservation += parts[parts.length - 1];
        }
      } else if (line.startsWith('#=GC cons')) {
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
          
          if (!seqMap[seqName]) {
            seqMap[seqName] = '';
          }
          
          seqMap[seqName] += seqSegment;
        }
      }
    }
    
    // Convert to array
    for (const [name, sequence] of Object.entries(seqMap)) {
      sequences.push({
        name,
        sequence,
        length: sequence.length
      });
    }
    
    return { sequences, conservation };
  };
  
  // Get color for residue based on property
  const getResidueColor = (residue, scheme) => {
    const residue_upper = residue.toUpperCase();
    
    // Default color
    if (!residue_upper || residue_upper === '-' || residue_upper === '.') {
      return '#cccccc';
    }
    
    // Different color schemes
    const schemes = {
      default: {
        A: '#77dd88', C: '#99ee66', D: '#55bb33', E: '#55bb33',
        F: '#9999ff', G: '#77dd88', H: '#7070ff', I: '#66bbff',
        K: '#ffcc77', L: '#66bbff', M: '#66bbff', N: '#55bb33',
        P: '#eeaaaa', Q: '#55bb33', R: '#ffcc77', S: '#ff4455',
        T: '#ff4455', V: '#66bbff', W: '#9999ff', Y: '#9999ff'
      },
      hydrophobicity: {
        A: '#ad0052', C: '#0000ff', D: '#0000ff', E: '#0000ff',
        F: '#ad0052', G: '#0000ff', H: '#0000ff', I: '#ad0052',
        K: '#0000ff', L: '#ad0052', M: '#ad0052', N: '#0000ff',
        P: '#0000ff', Q: '#0000ff', R: '#0000ff', S: '#0000ff',
        T: '#0000ff', V: '#ad0052', W: '#ad0052', Y: '#0000ff'
      },
      charge: {
        A: '#000000', C: '#000000', D: '#ff0000', E: '#ff0000',
        F: '#000000', G: '#000000', H: '#0000ff', I: '#000000',
        K: '#0000ff', L: '#000000', M: '#000000', N: '#000000',
        P: '#000000', Q: '#000000', R: '#0000ff', S: '#000000',
        T: '#000000', V: '#000000', W: '#000000', Y: '#000000'
      },
      nucleotide: {
        A: '#a0a0ff', T: '#a0ffa0', G: '#ff7070', C: '#ffff70',
        U: '#a0ffa0', N: '#ffffff'
      },
      clustal: {
        // Clustal color scheme
        A: '#80a0f0', R: '#f01505', N: '#00ff00', D: '#c048c0',
        C: '#f08080', Q: '#00ff00', E: '#c048c0', G: '#f09048',
        H: '#15a4a4', I: '#80a0f0', L: '#80a0f0', K: '#f01505',
        M: '#80a0f0', F: '#80a0f0', P: '#ffff00', S: '#00ff00',
        T: '#00ff00', W: '#80a0f0', Y: '#15a4a4', V: '#80a0f0',
        B: '#fff', Z: '#fff', X: '#fff', '-': '#fff'
      },
      taylor: {
        A: '#ccff00', R: '#0000ff', N: '#cc00ff', D: '#ff0000',
        C: '#ffff00', Q: '#ff00cc', E: '#ff0066', G: '#ff9900',
        H: '#0066ff', I: '#66ff00', L: '#33ff00', K: '#6600ff',
        M: '#00ff00', F: '#00ff66', P: '#ffcc00', S: '#ff3300',
        T: '#ff6600', W: '#00ccff', Y: '#00ffcc', V: '#99ff00',
        B: '#fff', Z: '#fff', X: '#fff', '-': '#fff'
      }
    };
    
    return schemes[scheme][residue_upper] || '#cccccc';
  };
  
  // Get color for conservation symbols
  const getConservationColor = (symbol) => {
    switch(symbol) {
      case '*': return '#ff0000'; // Fully conserved
      case ':': return '#ffa500'; // Strongly similar
      case '.': return '#0000ff'; // Weakly similar
      default: return 'transparent';
    }
  };
  
  // Navigate sequence view
  const navigateTo = (start) => {
    const width = viewportEnd - viewportStart;
    const newStart = Math.max(0, start);
    const newEnd = Math.min(sequences[0].length, newStart + width);
    
    setViewportStart(newStart);
    setViewportEnd(newEnd);
    setShowingFullWidth(newEnd === sequences[0].length);
  };
  
  // Adjust viewport width
  const adjustViewport = (delta) => {
    const maxLen = sequences[0]?.length || 100;
    const newWidth = Math.min(maxLen, Math.max(20, (viewportEnd - viewportStart) + delta));
    
    setViewportEnd(Math.min(maxLen, viewportStart + newWidth));
    setShowingFullWidth(viewportStart + newWidth >= maxLen);
  };
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div className="flex items-center">
          <label className="text-sm text-gray-600 mr-2">Color Scheme:</label>
          <select
            className="text-sm border rounded p-1"
            value={colorScheme}
            onChange={(e) => setColorScheme(e.target.value)}
          >
            <option value="clustal">Clustal</option>
            <option value="taylor">Taylor</option>
            <option value="hydrophobicity">Hydrophobicity</option>
            <option value="charge">Charge</option>
            <option value="nucleotide">Nucleotide</option>
          </select>
        </div>
        
        <div className="flex items-center space-x-2">
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm"
            onClick={() => navigateTo(viewportStart - Math.floor((viewportEnd - viewportStart) / 2))}
            disabled={viewportStart === 0}
          >
            ← Previous
          </button>
          
          <div className="flex items-center">
            <button
              className="p-1 text-gray-500 hover:bg-gray-100 rounded"
              onClick={() => setZoomLevel(Math.max(0.5, zoomLevel - 0.1))}
              title="Zoom Out"
            >
              <ZoomOut className="w-4 h-4" />
            </button>
            <div className="text-sm text-gray-600 w-16 text-center">{Math.round(zoomLevel * 100)}%</div>
            <button
              className="p-1 text-gray-500 hover:bg-gray-100 rounded"
              onClick={() => setZoomLevel(Math.min(3, zoomLevel + 0.1))}
              title="Zoom In"
            >
              <ZoomIn className="w-4 h-4" />
            </button>
          </div>
          
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm"
            onClick={() => navigateTo(viewportStart + Math.floor((viewportEnd - viewportStart) / 2))}
            disabled={showingFullWidth}
          >
            Next →
          </button>
          
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm"
            onClick={() => adjustViewport(20)}
            title="Show more columns"
          >
            +
          </button>
          
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm"
            onClick={() => adjustViewport(-20)}
            title="Show fewer columns"
            disabled={(viewportEnd - viewportStart) <= 20}
          >
            -
          </button>
        </div>
      </div>
      
      {/* Stats */}
      <div className="border-b p-2 text-sm text-gray-600">
        <span>Sequences: {sequences.length}</span>
        {sequences.length > 0 && (
          <>
            <span className="mx-2">|</span>
            <span>Alignment Length: {sequences[0].length} positions</span>
            <span className="mx-2">|</span>
            <span>
              Showing: {viewportStart + 1} - {viewportEnd}
              {!showingFullWidth && ` of ${sequences[0].length}`}
            </span>
          </>
        )}
      </div>
      
      {/* Alignment view */}
      <div className="flex-1 overflow-auto p-3">
        {sequences.length > 0 ? (
          <div className="font-mono" style={{ transform: `scale(${zoomLevel})`, transformOrigin: 'top left' }}>
            {/* Position ruler */}
            {settings.showRuler && (
              <div className="flex mb-1">
                <div className="inline-block w-32"></div>
                <div className="flex">
                  {Array.from({ length: Math.ceil((viewportEnd - viewportStart) / 10) }).map((_, i) => {
                    const pos = viewportStart + (i * 10);
                    return (
                      <span key={i} className="inline-block w-10 text-center text-xs text-gray-500">
                        {pos + 1}
                      </span>
                    );
                  })}
                </div>
              </div>
            )}
            
            {/* Sequences */}
            {sequences.map((seq, seqIndex) => (
              <div key={seqIndex} className="flex mb-1 hover:bg-gray-100">
                <div className="inline-block w-32 overflow-hidden text-ellipsis whitespace-nowrap text-sm text-right pr-2">
                  {seq.name}
                </div>
                <div className="flex">
                  {Array.from(seq.sequence.substring(viewportStart, viewportEnd)).map((residue, i) => (
                    <span 
                      key={i} 
                      className="inline-block text-center"
                      style={{ 
                        backgroundColor: getResidueColor(residue, colorScheme),
                        width: `${settings.fontSize * 0.8}px`,
                        height: `${settings.fontSize * 1.2}px`,
                        textAlign: 'center',
                        color: colorScheme === 'hydrophobicity' || colorScheme === 'charge' ? 'white' : 'black',
                        lineHeight: `${settings.fontSize * 1.2}px`
                      }}
                      title={`${seq.name}: ${viewportStart + i + 1}`}
                    >
                      {residue}
                    </span>
                  ))}
                </div>
              </div>
            ))}
            
            {/* Conservation line */}
            {conservation.length > 0 && (
              <div className="flex mt-1">
                <div className="inline-block w-32 text-sm text-right pr-2">Conservation</div>
                <div className="flex">
                  {Array.from(conservation.substring(viewportStart, viewportEnd)).map((symbol, i) => (
                    <span 
                      key={i} 
                      className="inline-block text-center"
                      style={{ 
                        color: getConservationColor(symbol),
                        fontWeight: 'bold',
                        width: `${settings.fontSize * 0.8}px`,
                        height: `${settings.fontSize * 1.2}px`,
                        textAlign: 'center',
                        lineHeight: `${settings.fontSize * 1.2}px`
                      }}
                    >
                      {symbol}
                    </span>
                  ))}
                </div>
              </div>
            )}
          </div>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No alignment data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Conservation viewer for multiple sequence alignments
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const ConservationViewer = ({ content, data, settings }) => {
  const chartRef = useRef(null);
  const [sequences, setSequences] = useState([]);
  const [conservation, setConservation] = useState([]);
  const [windowSize, setWindowSize] = useState(1);
  
  useEffect(() => {
    // Extract sequences based on format
    if (content) {
      let extractedSequences = [];
      
      if (data?.format === 'clustal') {
        const result = extractClustalData(content);
        extractedSequences = result.sequences;
      } else if (data?.format === 'stockholm') {
        const result = extractStockholmData(content);
        extractedSequences = result.sequences;
      }
      
      setSequences(extractedSequences);
      
      // Calculate conservation scores
      if (extractedSequences.length > 0) {
        const conservationScores = calculateConservation(extractedSequences);
        setConservation(conservationScores);
      }
    }
  }, [content, data]);
  
  // Extract data from CLUSTAL or Stockholm format
  const extractClustalData = (content) => {
    // Implementation similar to AlignmentViewer
    // ...
    const sequences = [];
    
    const lines = content.split('\n');
    
    // Find the header line
    const headerLineIndex = lines.findIndex(line => line.includes('CLUSTAL'));
    
    if (headerLineIndex !== -1) {
      const seqMap = {};
      
      // Parse sequences
      for (let i = headerLineIndex + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        
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
          name,
          sequence,
          length: sequence.length
        });
      }
    }
    
    return { sequences };
  };
  
  const extractStockholmData = (content) => {
    // Implementation similar to AlignmentViewer
    // ...
    const sequences = [];
    
    const lines = content.split('\n');
    
    const seqMap = {};
    
    for (const line of lines) {
      if (!line.startsWith('#') && !line.startsWith('//') && line.trim()) {
        // Sequence line
        const parts = line.trim().split(/\s+/);
        
        if (parts.length >= 2) {
          const seqName = parts[0];
          const seqSegment = parts[1];
          
          if (!seqMap[seqName]) {
            seqMap[seqName] = '';
          }
          
          seqMap[seqName] += seqSegment;
        }
      }
    }
    
    // Convert to array
    for (const [name, sequence] of Object.entries(seqMap)) {
      sequences.push({
        name,
        sequence,
        length: sequence.length
      });
    }
    
    return { sequences };
  };
  
  // Calculate conservation scores
  const calculateConservation = (sequences) => {
    if (sequences.length < 2) return [];
    
    const alignmentLength = sequences[0].length;
    const conservationScores = Array(alignmentLength).fill(0);
    
    // For each position in alignment
    for (let pos = 0; pos < alignmentLength; pos++) {
      // Count residues at this position
      const counts = {};
      let totalResidues = 0;
      
      // Count each residue
      for (const seq of sequences) {
        const residue = seq.sequence[pos];
        if (residue !== '-' && residue !== '.') {
          counts[residue] = (counts[residue] || 0) + 1;
          totalResidues++;
        }
      }
      
      // Skip positions with no residues
      if (totalResidues === 0) continue;
      
      // Calculate entropy-based conservation score
      let entropy = 0;
      for (const [_, count] of Object.entries(counts)) {
        const frequency = count / totalResidues;
        entropy -= frequency * Math.log(frequency);
      }
      
      // Normalize to 0-1 scale (1 is fully conserved)
      const maxEntropy = Math.log(Math.min(sequences.length, 20)); // 20 amino acids max
      const normalizedScore = 1 - (entropy / maxEntropy);
      
      conservationScores[pos] = normalizedScore;
    }
    
    return conservationScores;
  };
  
  // Apply smoothing window to conservation scores
  const smoothedConservation = useMemo(() => {
    if (conservation.length === 0 || windowSize <= 1) return conservation;
    
    const halfWindow = Math.floor(windowSize / 2);
    const smoothed = [];
    
    for (let i = 0; i < conservation.length; i++) {
      let sum = 0;
      let count = 0;
      
      for (let j = Math.max(0, i - halfWindow); j <= Math.min(conservation.length - 1, i + halfWindow); j++) {
        sum += conservation[j];
        count++;
      }
      
      smoothed.push(sum / count);
    }
    
    return smoothed;
  }, [conservation, windowSize]);
  
  // Render conservation plot
  useEffect(() => {
    if (!chartRef.current || conservation.length === 0) return;
    
    // Clear previous chart
    d3.select(chartRef.current).selectAll("*").remove();
    
    // Set up dimensions
    const margin = { top: 20, right: 30, bottom: 40, left: 50 };
    const width = chartRef.current.clientWidth - margin.left - margin.right;
    const height = chartRef.current.clientHeight - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(chartRef.current)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // X scale
    const x = d3.scaleLinear()
      .domain([1, conservation.length])
      .range([0, width]);
    
    // Y scale
    const y = d3.scaleLinear()
      .domain([0, 1])
      .range([height, 0]);
    
    // Draw axes
    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(x).ticks(10))
      .append("text")
        .attr("x", width / 2)
        .attr("y", margin.bottom - 5)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Position");
    
    svg.append("g")
      .call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".0%")))
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left + 15)
        .attr("x", -height / 2)
        .attr("fill", "black")
        .style("text-anchor", "middle")
        .text("Conservation");
    
    // Add conservation zones
    svg.append("rect")
      .attr("x", 0)
      .attr("y", y(1))
      .attr("width", width)
      .attr("height", y(0.8) - y(1))
      .attr("fill", "#c3e6c3")
      .attr("opacity", 0.3);
    
    svg.append("rect")
      .attr("x", 0)
      .attr("y", y(0.8))
      .attr("width", width)
      .attr("height", y(0.6) - y(0.8))
      .attr("fill", "#ffffbf")
      .attr("opacity", 0.3);
    
    svg.append("rect")
      .attr("x", 0)
      .attr("y", y(0.6))
      .attr("width", width)
      .attr("height", y(0) - y(0.6))
      .attr("fill", "#fdae61")
      .attr("opacity", 0.3);
    
    // Get data
    const data = smoothedConservation.map((score, i) => ({
      position: i + 1,
      score
    }));
    
    // Add the line
    svg.append("path")
      .datum(data)
      .attr("fill", "none")
      .attr("stroke", "#69b3a2")
      .attr("stroke-width", 2)
      .attr("d", d3.line()
        .x(d => x(d.position))
        .y(d => y(d.score))
      );
    
    // Add the area
    svg.append("path")
      .datum(data)
      .attr("fill", "#69b3a2")
      .attr("fill-opacity", 0.3)
      .attr("stroke", "none")
      .attr("d", d3.area()
        .x(d => x(d.position))
        .y0(y(0))
        .y1(d => y(d.score))
      );
    
    // Add tooltip
    const tooltip = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0)
      .style("background-color", "white")
      .style("border", "solid")
      .style("border-width", "1px")
      .style("border-radius", "5px")
      .style("padding", "5px")
      .style("position", "absolute");
    
    // Add interactive points
    svg.selectAll("circle")
      .data(data)
      .join("circle")
        .attr("cx", d => x(d.position))
        .attr("cy", d => y(d.score))
        .attr("r", 0) // Start with radius 0
        .attr("fill", d => {
          if (d.score >= 0.8) return "#1a9641";
          if (d.score >= 0.6) return "#a6611a";
          return "#d7191c";
        })
        .on("mouseover", function(d) {
          d3.select(this).transition()
            .duration(200)
            .attr("r", 5);
            
          tooltip.transition()
            .duration(200)
            .style("opacity", .9);
          
          // Get residues at this position
          let residues = '';
          if (sequences.length > 0) {
            residues = sequences.map(seq => 
              `${seq.name}: ${seq.sequence[d.position - 1]}`
            ).join('<br>');
          }
          
          tooltip.html(`Position: ${d.position}<br>Conservation: ${(d.score * 100).toFixed(1)}%<br><hr>${residues}`)
            .style("left", (event.pageX + 10) + "px")
            .style("top", (event.pageY - 28) + "px");
        })
        .on("mouseout", function(d) {
          d3.select(this).transition()
            .duration(200)
            .attr("r", 0);
            
          tooltip.transition()
            .duration(500)
            .style("opacity", 0);
        });
    
    // Add title
    svg.append("text")
      .attr("x", width / 2)
      .attr("y", -5)
      .attr("text-anchor", "middle")
      .style("font-size", "14px")
      .text(`Conservation Profile ${windowSize > 1 ? `(Window Size: ${windowSize})` : ''}`);
    
    // Return cleanup function
    return () => {
      d3.select("body").selectAll(".tooltip").remove();
    };
  }, [smoothedConservation, conservation, sequences, windowSize]);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div className="flex items-center">
          <label className="text-sm text-gray-600 mr-2">Smoothing Window:</label>
          <input
            type="range"
            min="1"
            max="15"
            step="2"
            value={windowSize}
            onChange={(e) => setWindowSize(parseInt(e.target.value))}
            className="w-32"
          />
          <span className="ml-2 text-sm">{windowSize}</span>
        </div>
        
        <div className="text-sm text-gray-600">
          <span>Sequences: {sequences.length}</span>
          {sequences.length > 0 && (
            <>
              <span className="mx-2">|</span>
              <span>Alignment Length: {sequences[0].length} positions</span>
            </>
          )}
        </div>
      </div>
      
      {/* Plot */}
      <div className="flex-1" ref={chartRef}>
        {conservation.length === 0 && (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No conservation data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Phylogenetic tree viewer component
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const PhylogeneticTreeViewer = ({ content, data, settings }) => {
  const svgRef = useRef(null);
  const [treeLayout, setTreeLayout] = useState('radial');
  const [showBranchLengths, setShowBranchLengths] = useState(true);
  const [showLabels, setShowLabels] = useState(true);
  
  useEffect(() => {
    if (!svgRef.current || !content) return;
    
    // Clear previous tree
    d3.select(svgRef.current).selectAll("*").remove();
    
    // Parse Newick format (simplified)
    const newickTree = extractNewickTree(content, data.format);
    
    if (!newickTree) {
      console.error("Failed to extract Newick tree");
      return;
    }
    
    // Set up dimensions
    const width = svgRef.current.clientWidth;
    const height = svgRef.current.clientHeight;
    const margin = 50;
    
    const svg = d3.select(svgRef.current)
      .append("svg")
        .attr("width", width)
        .attr("height", height)
      .append("g")
        .attr("transform", `translate(${margin}, ${margin})`);
    
    // Parse the Newick tree
    const root = d3.hierarchy(parseNewick(newickTree));
    
    // Compute tree layout
    let treeFunc;
    if (treeLayout === 'radial') {
      treeFunc = d3.tree()
        .size([2 * Math.PI, Math.min(width, height) / 2 - margin * 2])
        .separation((a, b) => (a.parent === b.parent ? 1 : 2) / a.depth);
    } else {
      // Rectangular layout
      treeFunc = d3.tree()
        .size([height - margin * 2, width - margin * 2]);
    }
    
    treeFunc(root);
    
    // Draw the tree based on layout
    if (treeLayout === 'radial') {
      drawRadialTree(svg, root, width, height, showBranchLengths, showLabels);
    } else {
      drawRectangularTree(svg, root, width, height, showBranchLengths, showLabels);
    }
    
  }, [content, data, treeLayout, showBranchLengths, showLabels]);
  
  // Extract Newick tree from different formats
  const extractNewickTree = (content, format) => {
    if (format === 'newick') {
      return content.trim();
    } else if (format === 'nexus') {
      // Extract tree from NEXUS format
      const treeMatch = content.match(/TREE\s+\w+\s*=\s*([^;]+;)/i);
      if (treeMatch) {
        return treeMatch[1].trim();
      }
    }
    return null;
  };
  
  // Parse Newick string (simplified)
  const parseNewick = (newick) => {
    // Simple recursive parsing of Newick format
    // This is a simplified version; a real implementation would be more complex
    
    // Remove whitespace
    newick = newick.replace(/\s+/g, '');
    
    // Function to parse a node
    const parseNode = (str, startPos) => {
      let pos = startPos || 0;
      const node = { name: '', children: [], length: null };
      
      // Check if we're at a leaf or internal node
      if (str[pos] === '(') {
        // Internal node - parse children
        pos++;
        let done = false;
        
        while (!done) {
          const childResult = parseNode(str, pos);
          pos = childResult.pos;
          node.children.push(childResult.node);
          
          if (str[pos] === ',') {
            pos++; // Move past comma
          } else if (str[pos] === ')') {
            pos++; // Move past closing paren
            done = true;
          }
        }
        
        // Check for node name
        let nameStr = '';
        while (pos < str.length && str[pos] !== ':' && str[pos] !== ',' && str[pos] !== ';' && str[pos] !== ')') {
          nameStr += str[pos];
          pos++;
        }
        node.name = nameStr;
      } else {
        // Leaf node - parse name
        let nameStr = '';
        while (pos < str.length && str[pos] !== ':' && str[pos] !== ',' && str[pos] !== ';' && str[pos] !== ')') {
          nameStr += str[pos];
          pos++;
        }
        node.name = nameStr;
      }
      
      // Check for branch length
      if (pos < str.length && str[pos] === ':') {
        pos++; // Move past colon
        let lengthStr = '';
        while (pos < str.length && str[pos] !== ',' && str[pos] !== ';' && str[pos] !== ')') {
          lengthStr += str[pos];
          pos++;
        }
        node.length = parseFloat(lengthStr);
      }
      
      return { node, pos };
    };
    
    // Parse the tree, removing trailing semicolon
    const result = parseNode(newick.endsWith(';') ? newick.slice(0, -1) : newick);
    return result.node;
  };
  
  // Draw radial tree
  const drawRadialTree = (svg, root, width, height, showBranchLengths, showLabels) => {
    // Draw links
    svg.append("g")
      .attr("fill", "none")
      .attr("stroke", "#555")
      .attr("stroke-opacity", 0.4)
      .attr("stroke-width", 1.5)
      .selectAll("path")
      .data(root.links())
      .join("path")
      .attr("d", d3.linkRadial()
        .angle(d => d.x)
        .radius(d => d.y)
      );
    
    // Draw nodes
    const nodes = svg.append("g")
      .selectAll("g")
      .data(root.descendants())
      .join("g")
      .attr("transform", d => `
        translate(${d.y * Math.sin(d.x)},${d.y * Math.cos(d.x) * -1})
      `);
    
    nodes.append("circle")
      .attr("fill", d => d.children ? "#555" : "#999")
      .attr("r", 3);
    
    // Add labels
    if (showLabels) {
      nodes.append("text")
        .attr("dy", "0.31em")
        .attr("x", d => d.x < Math.PI ? 6 : -6)
        .attr("text-anchor", d => d.x < Math.PI ? "start" : "end")
        .attr("transform", d => d.x >= Math.PI ? "rotate(180)" : null)
        .text(d => d.data.name)
        .clone(true).lower()
        .attr("stroke", "white");
    }
    
    // Add branch lengths
    if (showBranchLengths) {
      nodes.filter(d => d.data.length !== null && d.data.length !== undefined && d.parent)
        .append("text")
        .attr("dy", "-0.31em")
        .attr("x", d => d.x < Math.PI ? 6 : -6)
        .attr("text-anchor", d => d.x < Math.PI ? "start" : "end")
        .attr("transform", d => d.x >= Math.PI ? "rotate(180)" : null)
        .text(d => d.data.length.toFixed(3))
        .attr("font-size", "8px")
        .attr("fill", "#777");
    }
  };
  
  // Draw rectangular tree
  const drawRectangularTree = (svg, root, width, height, showBranchLengths, showLabels) => {
    // Draw links
    svg.append("g")
      .attr("fill", "none")
      .attr("stroke", "#555")
      .attr("stroke-opacity", 0.4)
      .attr("stroke-width", 1.5)
      .selectAll("path")
      .data(root.links())
      .join("path")
      .attr("d", d3.linkHorizontal()
        .x(d => d.y)
        .y(d => d.x)
      );
    
    // Draw nodes
    const nodes = svg.append("g")
      .selectAll("g")
      .data(root.descendants())
      .join("g")
      .attr("transform", d => `translate(${d.y},${d.x})`);
    
    nodes.append("circle")
      .attr("fill", d => d.children ? "#555" : "#999")
      .attr("r", 3);
    
    // Add labels
    if (showLabels) {
      nodes.append("text")
        .attr("dy", "0.31em")
        .attr("x", d => d.children ? -6 : 6)
        .attr("text-anchor", d => d.children ? "end" : "start")
        .text(d => d.data.name)
        .clone(true).lower()
        .attr("stroke", "white");
    }
    
    // Add branch lengths
    if (showBranchLengths) {
      svg.append("g")
        .selectAll("text")
        .data(root.links())
        .join("text")
        .attr("x", d => (d.source.y + d.target.y) / 2)
        .attr("y", d => (d.source.x + d.target.x) / 2 - 5)
        .attr("text-anchor", "middle")
        .text(d => d.target.data.length !== null ? d.target.data.length.toFixed(3) : "")
        .attr("font-size", "8px")
        .attr("fill", "#777");
    }
  };
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div className="flex items-center">
          <label className="text-sm text-gray-600 mr-2">Layout:</label>
          <select
            className="text-sm border rounded p-1"
            value={treeLayout}
            onChange={(e) => setTreeLayout(e.target.value)}
          >
            <option value="radial">Radial</option>
            <option value="rectangular">Rectangular</option>
          </select>
        </div>
        
        <div className="flex items-center space-x-4">
          <label className="flex items-center">
            <input
              type="checkbox"
              checked={showBranchLengths}
              onChange={() => setShowBranchLengths(!showBranchLengths)}
              className="mr-1"
            />
            <span className="text-sm text-gray-600">Branch Lengths</span>
          </label>
          
          <label className="flex items-center">
            <input
              type="checkbox"
              checked={showLabels}
              onChange={() => setShowLabels(!showLabels)}
              className="mr-1"
            />
            <span className="text-sm text-gray-600">Labels</span>
          </label>
        </div>
      </div>
      
      {/* Tree visualization */}
      <div className="flex-1" ref={svgRef}>
        {!content && (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No phylogenetic tree data available</p>
          </div>
        )}
      </div>
    </div>
  );
};

/**
 * Simple text viewer component
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
const TextViewer = ({ content, settings }) => {
  const [fontSize, setFontSize] = useState(settings.fontSize || 12);
  
  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="border-b p-2 flex justify-between items-center">
        <div className="flex items-center">
          <label className="text-sm text-gray-600 mr-2">Font Size:</label>
          <input
            type="range"
            min="8"
            max="18"
            step="1"
            value={fontSize}
            onChange={(e) => setFontSize(parseInt(e.target.value))}
            className="w-32"
          />
          <span className="ml-2 text-sm">{fontSize}px</span>
        </div>
        
        <div className="flex items-center">
          <button
            className="px-2 py-1 bg-gray-200 rounded text-sm flex items-center"
            onClick={() => {
              // Create download
              const blob = new Blob([content], { type: 'text/plain' });
              const url = URL.createObjectURL(blob);
              const a = document.createElement('a');
              a.href = url;
              a.download = 'file.txt';
              document.body.appendChild(a);
              a.click();
              document.body.removeChild(a);
              URL.revokeObjectURL(url);
            }}
          >
            <Download className="w-4 h-4 mr-1" />
            Download
          </button>
        </div>
      </div>
      
      {/* Text content */}
      <div className="flex-1 overflow-auto p-4">
        {content ? (
          <pre className="whitespace-pre-wrap" style={{ fontSize: `${fontSize}px` }}>
            {content}
          </pre>
        ) : (
          <div className="flex items-center justify-center h-full">
            <p className="text-gray-500">No text content available</p>
          </div>
        )}
      </div>
    </div>
  );
};

export default VisualizationComponent;
