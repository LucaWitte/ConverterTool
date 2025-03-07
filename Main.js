import React, { useState, useEffect, useRef } from 'react';
import { AlertCircle, Download, FileText, Upload, Check, ChevronDown, RefreshCw, FileCode, List, Settings, HelpCircle, Globe, Search, Database, Bookmark, Clock, BarChart2, Scissors, Zap, Layers, Save, Star } from 'lucide-react';

// Import modules
import { FormatDefinitions, FormatCategories, ConversionMatrix } from './FormatDefinitions';
import { DatabaseResources, DatabaseUrlToFormat, ExampleFiles } from './DatabaseResources';
import { detectFileFormat, validateFormat } from './FormatDetection';
import { parseFasta, parseFastq, parseGenBank, parsePDB, parseClustal } from './FormatParsers';
import { convertFile, performConversion, formatSequenceForFasta } from './FormatConverters';
import { performAnalysis, AnalysisTools } from './AnalysisTools';
import { FileDropZone } from './UIComponents';
import { Notification } from './UIComponents';
import { TabNavigation } from './UIComponents';
import { UrlLoader } from './UrlLoader';
import { DatabaseSearch } from './DatabaseSearch';
import { ConversionOptions } from './ConversionOptions';

const BioFormatConverter = () => {
  // State management
  const [uploadedFile, setUploadedFile] = useState(null);
  const [fileContent, setFileContent] = useState('');
  const [fileType, setFileType] = useState('');
  const [detectedFormat, setDetectedFormat] = useState('');
  const [targetFormat, setTargetFormat] = useState('');
  const [previewData, setPreviewData] = useState(null);
  const [convertedData, setConvertedData] = useState(null);
  const [isConverting, setIsConverting] = useState(false);
  const [error, setError] = useState(null);
  const [errorDetails, setErrorDetails] = useState(null); // Added for detailed error tracking
  const [showFormatInfo, setShowFormatInfo] = useState(false);
  const [showAdvancedOptions, setShowAdvancedOptions] = useState(false);
  const [activeTab, setActiveTab] = useState('upload');
  const [urlInput, setUrlInput] = useState('');
  const [isLoadingUrl, setIsLoadingUrl] = useState(false);
  const [databaseSearch, setDatabaseSearch] = useState('');
  const [searchResults, setSearchResults] = useState([]);
  const [isSearching, setIsSearching] = useState(false);
  const [savedJobs, setSavedJobs] = useState([]);
  const [analysisResults, setAnalysisResults] = useState(null);
  const [showAnalysisOptions, setShowAnalysisOptions] = useState(false);
  const [conversionOptions, setConversionOptions] = useState({
    preserveHeaders: true,
    formatSequence: true,
    validateOutput: true,
    lineLength: 60,
    includeQuality: true,
    addSourceComment: true,
    normalizeSequence: false,
    translateToProtein: false,
    removeGaps: false,
    includeFeatures: true
  });
  const [processingProgress, setProcessingProgress] = useState(0);
  const [showHelpPanel, setShowHelpPanel] = useState(false);
  const [selectedExample, setSelectedExample] = useState(null);
  const [dragActive, setDragActive] = useState(false);
  const [savedFormats, setSavedFormats] = useState([]);
  const [showNotification, setShowNotification] = useState(false);
  const [notificationMessage, setNotificationMessage] = useState('');
  const [jobHistory, setJobHistory] = useState([]);
  
  const fileInputRef = useRef(null);
  
  // Initialize effect for browser storage
  useEffect(() => {
    try {
      // Load saved jobs from localStorage if available
      const savedJobsData = localStorage.getItem('bioConverter_savedJobs');
      if (savedJobsData) {
        setSavedJobs(JSON.parse(savedJobsData));
      }
      
      // Load job history
      const historyData = localStorage.getItem('bioConverter_history');
      if (historyData) {
        setJobHistory(JSON.parse(historyData));
      }
      
      // Load saved format preferences
      const savedFormatsData = localStorage.getItem('bioConverter_savedFormats');
      if (savedFormatsData) {
        setSavedFormats(JSON.parse(savedFormatsData));
      }
    } catch (e) {
      console.error("Error loading data from localStorage:", e);
      setError("Failed to load saved data. Local storage might be corrupted or unavailable.");
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "localStorage initialization"
      });
    }
  }, []);
  
  // Save jobs to localStorage when they change
  useEffect(() => {
    if (savedJobs.length > 0) {
      try {
        localStorage.setItem('bioConverter_savedJobs', JSON.stringify(savedJobs));
      } catch (e) {
        console.error("Error saving jobs to localStorage:", e);
        showNotificationMessage("Failed to save jobs to local storage. Your browser storage might be full or restricted.");
      }
    }
  }, [savedJobs]);
  
  // Save history to localStorage when it changes
  useEffect(() => {
    if (jobHistory.length > 0) {
      try {
        localStorage.setItem('bioConverter_history', JSON.stringify(jobHistory));
      } catch (e) {
        console.error("Error saving history to localStorage:", e);
        showNotificationMessage("Failed to save history to local storage. Your browser storage might be full or restricted.");
      }
    }
  }, [jobHistory]);
  
  // Reset application state
  const resetState = () => {
    setUploadedFile(null);
    setFileContent('');
    setFileType('');
    setDetectedFormat('');
    setTargetFormat('');
    setPreviewData(null);
    setConvertedData(null);
    setError(null);
    setErrorDetails(null);
    setProcessingProgress(0);
    setUrlInput('');
    setAnalysisResults(null);
  };
  
  // Handle file upload
  const handleFileUpload = (event) => {
    try {
      const file = event.target.files[0];
      if (!file) {
        setError("No file selected");
        return;
      }
      
      // Check file size
      if (file.size > 10 * 1024 * 1024) { // 10MB limit
        setError("File is too large. Maximum file size is 10MB.");
        return;
      }
      
      // Check file extension for basic validation
      const fileExtension = file.name.split('.').pop().toLowerCase();
      const knownExtensions = [].concat(...Object.values(FormatDefinitions).map(format => 
        format.extensions.map(ext => ext.substring(1).toLowerCase())
      ));
      
      if (!knownExtensions.includes(fileExtension)) {
        showNotificationMessage(`Warning: File extension .${fileExtension} is not recognized. Will attempt to detect format from content.`);
      }
      
      processUploadedFile(file);
    } catch (e) {
      console.error("File upload error:", e);
      setError(`Failed to upload file: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "handleFileUpload"
      });
      setProcessingProgress(0);
    }
  };
  
  // Process the uploaded file
  const processUploadedFile = (file) => {
    resetState();
    setUploadedFile(file);
    
    // Start reading and processing the file
    const reader = new FileReader();
    
    reader.onload = (e) => {
      try {
        // Simulate progressive loading for large files
        setProcessingProgress(50);
        
        const content = e.target.result;
        if (!content || content.trim() === '') {
          throw new Error("File is empty");
        }
        
        setFileContent(content);
        
        // Add a small delay to simulate processing for large files
        setTimeout(() => {
          try {
            // Auto-detect format based on file content
            const detectedFormat = detectFileFormat(content, file.name, FormatDefinitions);
            setFileType(detectedFormat);
            setDetectedFormat(detectedFormat);
            
            setProcessingProgress(75);
            
            setTimeout(() => {
              try {
                if (detectedFormat) {
                  generatePreview(content, detectedFormat);
                  setProcessingProgress(100);
                  
                  // Add to job history
                  addToHistory({
                    type: 'upload',
                    format: detectedFormat,
                    fileName: file.name,
                    timestamp: new Date().toISOString()
                  });
                } else {
                  setError("Unable to detect file format. Please select manually.");
                  setProcessingProgress(100);
                }
              } catch (e) {
                console.error("Error in preview generation:", e);
                setError(`Failed to process file: ${e.message}`);
                setErrorDetails({
                  message: e.message,
                  stack: e.stack,
                  source: "processUploadedFile - preview generation"
                });
                setProcessingProgress(0);
              }
            }, 300);
          } catch (e) {
            console.error("Error in format detection:", e);
            setError(`Failed to detect file format: ${e.message}`);
            setErrorDetails({
              message: e.message,
              stack: e.stack,
              source: "processUploadedFile - format detection"
            });
            setProcessingProgress(0);
          }
        }, 200);
      } catch (e) {
        console.error("Error processing file content:", e);
        setError(`Failed to process file content: ${e.message}`);
        setErrorDetails({
          message: e.message,
          stack: e.stack,
          source: "processUploadedFile - content processing"
        });
        setProcessingProgress(0);
      }
    };
    
    reader.onerror = (e) => {
      console.error("FileReader error:", e);
      setError("Error reading file. The file might be corrupted or inaccessible.");
      setErrorDetails({
        message: e.target.error?.message || "Unknown FileReader error",
        code: e.target.error?.code,
        source: "FileReader.onerror"
      });
      setProcessingProgress(0);
    };
    
    reader.onabort = () => {
      console.warn("File reading aborted");
      setError("File reading was aborted.");
      setProcessingProgress(0);
    };
    
    try {
      reader.readAsText(file);
    } catch (e) {
      console.error("Error initiating file read:", e);
      setError(`Failed to read file: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "reader.readAsText"
      });
      setProcessingProgress(0);
    }
  };
  
  // Generate preview data based on file format
  const generatePreview = (content, format) => {
    try {
      let preview;
      
      switch (format) {
        case "fasta":
          preview = parseFasta(content);
          break;
        case "fastq":
          preview = parseFastq(content);
          break;
        case "genbank":
          preview = parseGenBank(content);
          break;
        case "pdb":
          preview = parsePDB(content);
          break;
        case "clustal":
          preview = parseClustal(content);
          break;
        default:
          // Validate the format is supported
          if (!FormatDefinitions[format]) {
            throw new Error(`Format '${format}' is not supported for preview generation.`);
          }
          preview = { 
            format,
            error: "Preview not available for this format",
            message: `The format ${FormatDefinitions[format].name} is supported for conversion, but detailed preview is not available.`
          };
      }
      
      // Validate preview results
      if (!preview || (preview.error && !preview.message)) {
        throw new Error("Preview generation failed with an unknown error");
      }
      
      setPreviewData(preview);
      setError(null);
      setErrorDetails(null);
    } catch (err) {
      console.error("Preview generation error:", err);
      setError(`Error parsing file: ${err.message}`);
      setErrorDetails({
        message: err.message,
        stack: err.stack,
        source: "generatePreview",
        format: format
      });
      setPreviewData(null);
    }
  };
  
  // Load an example file
  const loadExampleFile = (formatKey) => {
    try {
      if (!ExampleFiles[formatKey]) {
        throw new Error(`Example file for format '${formatKey}' is not available.`);
      }
      
      resetState();
      setSelectedExample(formatKey);
      
      const example = ExampleFiles[formatKey];
      const file = new File([example.content], example.name, { type: 'text/plain' });
      
      setUploadedFile(file);
      setFileContent(example.content);
      setFileType(formatKey);
      setDetectedFormat(formatKey);
      
      generatePreview(example.content, formatKey);
      
      // Add to job history
      addToHistory({
        type: 'example',
        format: formatKey,
        fileName: example.name,
        timestamp: new Date().toISOString()
      });
      
      showNotificationMessage(`Loaded ${example.name} example file`);
    } catch (e) {
      console.error("Error loading example file:", e);
      setError(`Failed to load example file: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "loadExampleFile",
        formatKey: formatKey
      });
    }
  };
  
  // Handle drag events
  const handleDrag = (e) => {
    e.preventDefault();
    e.stopPropagation();
    
    if (e.type === "dragenter" || e.type === "dragover") {
      setDragActive(true);
    } else if (e.type === "dragleave") {
      setDragActive(false);
    }
  };
  
  // Handle drop event
  const handleDrop = (e) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
    
    try {
      if (e.dataTransfer.files && e.dataTransfer.files[0]) {
        const file = e.dataTransfer.files[0];
        
        // Check file size
        if (file.size > 10 * 1024 * 1024) { // 10MB limit
          setError("File is too large. Maximum file size is 10MB.");
          return;
        }
        
        processUploadedFile(file);
      } else {
        throw new Error("No valid files found in drop event");
      }
    } catch (err) {
      console.error("File drop error:", err);
      setError(`Failed to process dropped file: ${err.message}`);
      setErrorDetails({
        message: err.message,
        stack: err.stack,
        source: "handleDrop"
      });
    }
  };
  
  // Add job to history
  const addToHistory = (job) => {
    try {
      if (!job || !job.type) {
        throw new Error("Invalid job data for history");
      }
      
      setJobHistory(prev => {
        const newHistory = [
          { ...job, id: `job_${Date.now()}` },
          ...prev.slice(0, 19) // Keep only the 20 most recent jobs
        ];
        return newHistory;
      });
    } catch (e) {
      console.error("Error adding to history:", e);
      // Don't set error state as this is a non-critical operation
      // Just log the error and continue
    }
  };
  
  // Save current job
  const saveCurrentJob = () => {
    try {
      if (!uploadedFile || !fileType) {
        setError("No data to save.");
        return;
      }
      
      if (!fileContent) {
        throw new Error("File content is empty");
      }
      
      const jobToSave = {
        id: `saved_${Date.now()}`,
        name: uploadedFile.name,
        format: fileType,
        formatName: FormatDefinitions[fileType]?.name || "Unknown format",
        content: fileContent,
        savedAt: new Date().toISOString()
      };
      
      setSavedJobs(prev => [jobToSave, ...prev]);
      showNotificationMessage(`Saved ${uploadedFile.name}`);
    } catch (e) {
      console.error("Error saving job:", e);
      setError(`Failed to save job: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "saveCurrentJob"
      });
    }
  };
  
  // Load saved job
  const loadSavedJob = (job) => {
    try {
      if (!job || !job.content || !job.format) {
        throw new Error("Invalid saved job data");
      }
      
      resetState();
      
      setFileContent(job.content);
      setFileType(job.format);
      setDetectedFormat(job.format);
      
      // Validate the content matches the expected format
      const validationResult = validateFormat(job.content, job.format);
      if (!validationResult.valid) {
        showNotificationMessage(`Warning: Saved job may have format issues: ${validationResult.errors[0]}`);
      }
      
      generatePreview(job.content, job.format);
      
      // Create a file object
      const file = new File([job.content], job.name, { type: 'text/plain' });
      setUploadedFile(file);
      
      showNotificationMessage(`Loaded saved job: ${job.name}`);
    } catch (e) {
      console.error("Error loading saved job:", e);
      setError(`Failed to load saved job: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "loadSavedJob",
        jobId: job?.id
      });
    }
  };
  
  // Delete saved job
  const deleteSavedJob = (jobId) => {
    try {
      if (!jobId) {
        throw new Error("Invalid job ID for deletion");
      }
      
      setSavedJobs(prev => prev.filter(job => job.id !== jobId));
      showNotificationMessage("Saved job deleted");
    } catch (e) {
      console.error("Error deleting saved job:", e);
      showNotificationMessage(`Error deleting job: ${e.message}`);
    }
  };
  
  // Show notification
  const showNotificationMessage = (message) => {
    if (!message) return;
    
    setNotificationMessage(message);
    setShowNotification(true);
    
    // Hide after 3 seconds
    setTimeout(() => {
      setShowNotification(false);
    }, 3000);
  };
  
  // Handle conversion between formats
  const handleConversion = () => {
    try {
      if (!fileContent) {
        throw new Error("No file content to convert");
      }
      
      if (!fileType) {
        throw new Error("Source format not detected or selected");
      }
      
      if (!targetFormat) {
        throw new Error("Please select a target format for conversion");
      }
      
      // Validate source format is supported for conversion
      if (!ConversionMatrix[fileType]) {
        throw new Error(`Source format '${fileType}' does not support conversion`);
      }
      
      // Validate target format is supported for this source
      if (!ConversionMatrix[fileType].includes(targetFormat)) {
        throw new Error(`Conversion from ${fileType} to ${targetFormat} is not supported`);
      }
      
      setIsConverting(true);
      setError(null);
      setErrorDetails(null);
      setProcessingProgress(10);
      
      // Use the convertFile function from the FormatConverters module
      convertFile(
        fileContent, 
        fileType, 
        targetFormat, 
        conversionOptions,
        setProcessingProgress,
        (result) => {
          // Success callback
          if (!result || !result.content) {
            throw new Error("Conversion returned empty or invalid result");
          }
          
          setConvertedData(result);
          setIsConverting(false);
          setProcessingProgress(100);
          
          // Add to job history
          addToHistory({
            type: 'conversion',
            sourceFormat: fileType,
            targetFormat: targetFormat,
            fileName: uploadedFile ? uploadedFile.name : 'unknown',
            timestamp: new Date().toISOString()
          });
          
          showNotificationMessage(`Conversion complete: ${FormatDefinitions[fileType]?.name || fileType} → ${FormatDefinitions[targetFormat]?.name || targetFormat}`);
        },
        (errorMsg) => {
          // Error callback
          console.error("Conversion error:", errorMsg);
          setError(`Conversion failed: ${errorMsg}`);
          setErrorDetails({
            message: errorMsg,
            source: "convertFile error callback"
          });
          setIsConverting(false);
          setProcessingProgress(0);
        }
      );
    } catch (e) {
      console.error("Conversion initiation error:", e);
      setError(`Failed to start conversion: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "handleConversion"
      });
      setIsConverting(false);
      setProcessingProgress(0);
    }
  };
  
  // Handle analysis
  const handleAnalysis = (toolId) => {
    try {
      if (!fileContent || !fileType) {
        throw new Error("No data to analyze");
      }
      
      if (!toolId) {
        throw new Error("No analysis tool selected");
      }
      
      // Validate tool is supported for this format
      const availableTools = AnalysisTools[fileType] || [];
      if (!availableTools.some(tool => tool.id === toolId)) {
        throw new Error(`Analysis tool '${toolId}' is not supported for format '${fileType}'`);
      }
      
      setShowAnalysisOptions(false);
      setIsConverting(true);
      setProcessingProgress(10);
      
      // Use the performAnalysis function from the AnalysisTools module
      performAnalysis(
        toolId,
        fileContent,
        fileType,
        setProcessingProgress,
        (result) => {
          // Success callback
          if (!result) {
            throw new Error("Analysis returned empty or invalid result");
          }
          
          setAnalysisResults(result);
          setIsConverting(false);
          setProcessingProgress(100);
          
          // Add to job history
          addToHistory({
            type: 'analysis',
            tool: toolId,
            format: fileType,
            fileName: uploadedFile ? uploadedFile.name : 'unknown',
            timestamp: new Date().toISOString()
          });
          
          showNotificationMessage(`Analysis complete: ${result.title}`);
        },
        (errorMsg) => {
          // Error callback
          console.error("Analysis error:", errorMsg);
          setError(`Analysis failed: ${errorMsg}`);
          setErrorDetails({
            message: errorMsg,
            source: "performAnalysis error callback",
            toolId: toolId
          });
          setAnalysisResults(null);
          setIsConverting(false);
          setProcessingProgress(0);
        }
      );
    } catch (e) {
      console.error("Analysis initiation error:", e);
      setError(`Failed to start analysis: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "handleAnalysis",
        toolId: toolId
      });
      setIsConverting(false);
      setProcessingProgress(0);
    }
  };
  
  // Download converted file
  const downloadConvertedFile = () => {
    try {
      if (!convertedData || !convertedData.content) {
        throw new Error("No converted data available for download");
      }
      
      const formatKey = convertedData.format;
      if (!formatKey || !FormatDefinitions[formatKey]) {
        throw new Error("Invalid format for converted data");
      }
      
      const formatInfo = FormatDefinitions[formatKey];
      const extension = formatInfo.extensions[0];
      
      const blob = new Blob([convertedData.content], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      
      try {
        const a = document.createElement('a');
        a.href = url;
        a.download = `converted_file${extension}`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
      } finally {
        // Always revoke the object URL to avoid memory leaks
        URL.revokeObjectURL(url);
      }
    } catch (e) {
      console.error("Download error:", e);
      setError(`Failed to download file: ${e.message}`);
      setErrorDetails({
        message: e.message,
        stack: e.stack,
        source: "downloadConvertedFile"
      });
      showNotificationMessage(`Download failed: ${e.message}`);
    }
  };
  
  // Render format selection options
  const renderFormatOptions = () => {
    if (!fileType) return null;
    
    const compatibleFormats = ConversionMatrix[fileType] || [];
    
    return (
      <div className="mt-4">
        <h3 className="text-lg font-medium mb-2">Convert to:</h3>
        <select 
          className="w-full p-2 border rounded bg-white" 
          value={targetFormat}
          onChange={(e) => setTargetFormat(e.target.value)}
        >
          <option value="">Select target format</option>
          {compatibleFormats.map(format => (
            <option key={format} value={format}>
              {FormatDefinitions[format].name} ({FormatDefinitions[format].extensions.join(', ')})
            </option>
          ))}
        </select>
        
        {showAdvancedOptions && (
          <ConversionOptions 
            options={conversionOptions}
            setOptions={setConversionOptions}
            sourceFormat={fileType}
            targetFormat={targetFormat}
          />
        )}
        
        <div className="mt-2">
          <button
            className="text-sm text-blue-600 flex items-center"
            onClick={() => setShowAdvancedOptions(!showAdvancedOptions)}
          >
            <Settings className="w-4 h-4 mr-1" />
            {showAdvancedOptions ? 'Hide advanced options' : 'Show advanced options'}
          </button>
        </div>
        
        <button 
          className="mt-3 bg-blue-600 text-white p-2 rounded w-full flex justify-center items-center space-x-2 disabled:bg-blue-300"
          onClick={handleConversion}
          disabled={!targetFormat || isConverting}
        >
          {isConverting ? (
            <>
              <RefreshCw className="w-5 h-5 animate-spin" />
              <span>Converting... {processingProgress}%</span>
            </>
          ) : (
            <>
              <RefreshCw className="w-5 h-5" />
              <span>Convert</span>
            </>
          )}
        </button>
      </div>
    );
  };
  
  // Render the analysis tools available for the current format
  const renderAnalysisTools = () => {
    if (!fileType || !fileContent) return null;
    
    const tools = AnalysisTools[fileType] || [];
    
    if (tools.length === 0) return null;
    
    return (
      <div className="mt-4">
        <button
          className="flex items-center text-blue-600 text-sm"
          onClick={() => setShowAnalysisOptions(!showAnalysisOptions)}
        >
          <BarChart2 className="w-4 h-4 mr-1" />
          <span>{showAnalysisOptions ? 'Hide analysis tools' : 'Show analysis tools'}</span>
          <ChevronDown className={`w-4 h-4 ml-1 transform ${showAnalysisOptions ? 'rotate-180' : ''}`} />
        </button>
        
        {showAnalysisOptions && (
          <div className="mt-2 border rounded-md p-3 bg-gray-50">
            <h4 className="font-medium mb-2">Available Analysis Tools:</h4>
            <div className="space-y-2">
              {tools.map(tool => (
                <button
                  key={tool.id}
                  className="w-full text-left p-2 rounded hover:bg-blue-50 flex items-center"
                  onClick={() => handleAnalysis(tool.id)}
                  disabled={isConverting}
                >
                  <Zap className="w-4 h-4 mr-2 text-blue-500" />
                  <div>
                    <div className="font-medium">{tool.name}</div>
                    <div className="text-xs text-gray-600">{tool.description}</div>
                  </div>
                </button>
              ))}
            </div>
          </div>
        )}
      </div>
    );
  };
  
  // Render file preview based on format
  const renderPreview = () => {
    if (!previewData) return null;
    
    return (
      <div className="mt-6 border p-4 rounded bg-gray-50">
        <h3 className="text-lg font-medium mb-2">File Preview:</h3>
        
        {previewData.format === "fasta" && (
          <div>
            <p><strong>Format:</strong> FASTA</p>
            <p><strong>Sequences:</strong> {previewData.sequenceCount}</p>
            <p><strong>Total Length:</strong> {previewData.totalLength} bp</p>
            <p><strong>Type:</strong> {previewData.type === "nucleotide" ? "Nucleotide sequence" : "Protein sequence"}</p>
            
            {previewData.gcContent && (
              <p><strong>GC Content:</strong> {previewData.gcContent}</p>
            )}
            
            {previewData.sequences.length > 0 && (
              <div className="mt-2 border-t pt-2">
                <p><strong>First sequence:</strong></p>
                <p className="text-sm text-gray-600 truncate">{previewData.sequences[0].header}</p>
                <p className="text-sm font-mono whitespace-pre-wrap overflow-hidden" style={{maxHeight: '100px'}}>
                  {previewData.sequences[0].sequence}
                </p>
              </div>
            )}
          </div>
        )}
        
        {previewData.format === "fastq" && (
          <div>
            <p><strong>Format:</strong> FASTQ</p>
            <p><strong>Sequences:</strong> {previewData.sequenceCount}</p>
            <p><strong>Average Read Length:</strong> {previewData.avgReadLength} bp</p>
            <p><strong>Quality scores:</strong> min: {previewData.qualityStats.min}, max: {previewData.qualityStats.max}, avg: {previewData.qualityStats.avg}</p>
            
            {previewData.sequences.length > 0 && (
              <div className="mt-2 border-t pt-2">
                <p><strong>Sample read:</strong></p>
                <p className="text-sm text-gray-600 truncate">@{previewData.sequences[0].header}</p>
                <p className="text-sm font-mono">{previewData.sequences[0].sequence}</p>
                <p className="text-sm">+</p>
                <p className="text-sm font-mono">{previewData.sequences[0].quality}</p>
              </div>
            )}
          </div>
        )}
        
        {previewData.format === "genbank" && (
          <div>
            <p><strong>Format:</strong> GenBank</p>
            <p><strong>Sequence:</strong> {previewData.metadata.name || 'Unknown'}</p>
            <p><strong>Length:</strong> {previewData.sequenceLength} bp</p>
            <p><strong>Features:</strong> {previewData.featureCount}</p>
            
            {previewData.metadata.definition && (
              <p className="mt-2"><strong>Definition:</strong> {previewData.metadata.definition}</p>
            )}
            
            {previewData.metadata.organism && (
              <p><strong>Organism:</strong> {previewData.metadata.organism}</p>
            )}
            
            {previewData.featureTypes && Object.keys(previewData.featureTypes).length > 0 && (
              <div className="mt-2 border-t pt-2">
                <p><strong>Feature types:</strong></p>
                <ul className="list-disc list-inside text-sm">
                  {Object.entries(previewData.featureTypes).map(([type, count]) => (
                    <li key={type}>{type}: {count}</li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        )}
        
        {previewData.format === "pdb" && (
          <div>
            <p><strong>Format:</strong> PDB</p>
            {previewData.metadata.title && (
              <p><strong>Title:</strong> {previewData.metadata.title}</p>
            )}
            <p><strong>Atoms:</strong> {previewData.atomCount}</p>
            <p><strong>HETATM records:</strong> {previewData.hetatmCount}</p>
            <p><strong>Chains:</strong> {previewData.chainCount} ({previewData.chains.join(', ')})</p>
            
            {previewData.secondaryStructure && (
              <p><strong>Secondary structure:</strong> {previewData.secondaryStructure.helixCount} helices, {previewData.secondaryStructure.sheetCount} sheets</p>
            )}
            
            {previewData.residueTypes.length > 0 && (
              <div className="mt-2">
                <p><strong>Residue types:</strong> {previewData.residueTypes.slice(0, 10).join(', ')}{previewData.residueTypes.length > 10 ? '...' : ''}</p>
              </div>
            )}
            
            {previewData.atoms.length > 0 && (
              <div className="mt-2 border-t pt-2">
                <p><strong>Sample atoms:</strong></p>
                <div className="overflow-x-auto">
                  <table className="text-xs w-full mt-1">
                    <thead>
                      <tr className="border-b">
                        <th className="px-2 py-1 text-left">Atom</th>
                        <th className="px-2 py-1 text-left">Residue</th>
                        <th className="px-2 py-1 text-left">Chain</th>
                        <th className="px-2 py-1 text-left">Position</th>
                      </tr>
                    </thead>
                    <tbody>
                      {previewData.atoms.map((atom, index) => (
                        <tr key={index} className={index % 2 === 0 ? 'bg-gray-100' : ''}>
                          <td className="px-2 py-1">{atom.atomName}</td>
                          <td className="px-2 py-1">{atom.residueName} {atom.residueNum}</td>
                          <td className="px-2 py-1">{atom.chainId}</td>
                          <td className="px-2 py-1">({atom.x.toFixed(1)}, {atom.y.toFixed(1)}, {atom.z.toFixed(1)})</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}
          </div>
        )}
        
        {previewData.format === "clustal" && (
          <div>
            <p><strong>Format:</strong> Clustal</p>
            <p><strong>Sequences:</strong> {previewData.sequenceCount}</p>
            <p><strong>Alignment Length:</strong> {previewData.alignmentLength} positions</p>
            
            {previewData.identities && previewData.identities.length > 0 && (
              <p><strong>Identity:</strong> {previewData.identities[0].seq1} vs {previewData.identities[0].seq2}: {previewData.identities[0].identityPercent}</p>
            )}
            
            {previewData.conservedRegions && previewData.conservedRegions.length > 0 && (
              <div className="mt-2">
                <p><strong>Conserved Regions:</strong></p>
                <ul className="list-disc list-inside text-sm">
                  {previewData.conservedRegions.map((region, index) => (
                    <li key={index}>Position {region.start}-{region.end} (length: {region.length})</li>
                  ))}
                </ul>
              </div>
            )}
            
            {previewData.sequences.length > 0 && (
              <div className="mt-2 border-t pt-2">
                <p><strong>Alignment Preview:</strong></p>
                <div className="font-mono text-xs whitespace-pre-wrap overflow-x-auto">
                  {previewData.sequences.map((seq, index) => (
                    <div key={index} className="mt-1">
                      <span className="inline-block w-20 font-semibold">{seq.name}</span>
                      <span>{seq.sequence}</span>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}
      </div>
    );
  };
  
  // Render conversion result
  const renderConversionResult = () => {
    if (!convertedData) return null;
    
    return (
      <div className="mt-6 border p-4 rounded bg-green-50">
        <div className="flex justify-between items-start">
          <h3 className="text-lg font-medium mb-2 flex items-center">
            <Check className="w-5 h-5 text-green-600 mr-2" />
            Conversion Successful
          </h3>
          
          <button
            className="bg-green-600 text-white p-2 rounded flex items-center space-x-1"
            onClick={downloadConvertedFile}
          >
            <Download className="w-4 h-4" />
            <span>Download</span>
          </button>
        </div>
        
        <div className="mt-2">
          <p><strong>From:</strong> {FormatDefinitions[convertedData.originalFormat].name}</p>
          <p><strong>To:</strong> {FormatDefinitions[convertedData.format].name}</p>
          
          {convertedData.sequenceCount && (
            <p><strong>Sequences:</strong> {convertedData.sequenceCount}</p>
          )}
          
          {convertedData.totalLength && (
            <p><strong>Total Length:</strong> {convertedData.totalLength} bp</p>
          )}
          
          {convertedData.note && (
            <p className="mt-2 text-sm text-amber-700">{convertedData.note}</p>
          )}
        </div>
        
        <div className="mt-3 border-t pt-3">
          <h4 className="font-medium mb-1">Preview:</h4>
          <div className="bg-white border rounded p-2 max-h-52 overflow-auto">
            <pre className="text-xs">{convertedData.content.substring(0, 500)}{convertedData.content.length > 500 ? '...' : ''}</pre>
          </div>
        </div>
      </div>
    );
  };
  
  // Render analysis results
  const renderAnalysisResults = () => {
    if (!analysisResults) return null;
    
    return (
      <div className="mt-6 border p-4 rounded bg-blue-50">
        <h3 className="text-lg font-medium mb-2 flex items-center">
          <BarChart2 className="w-5 h-5 text-blue-600 mr-2" />
          {analysisResults.title}
        </h3>
        
        {/* Sequence Statistics Analysis */}
        {analysisResults.type === "sequence_stats" && (
          <div>
            <div className="grid grid-cols-2 gap-4 mb-4">
              <div className="bg-white p-3 rounded shadow-sm">
                <p className="font-medium mb-1">Overview</p>
                <p className="text-sm">Sequence Type: {analysisResults.data.sequenceType}</p>
                <p className="text-sm">Sequences: {analysisResults.data.stats.sequenceCount}</p>
                <p className="text-sm">Total Length: {analysisResults.data.stats.totalLength} bp</p>
                <p className="text-sm">Average Length: {analysisResults.data.stats.averageLength} bp</p>
              </div>
              
              <div className="bg-white p-3 rounded shadow-sm">
                <p className="font-medium mb-1">Length Statistics</p>
                <p className="text-sm">Min Length: {analysisResults.data.stats.minLength} bp</p>
                <p className="text-sm">Max Length: {analysisResults.data.stats.maxLength} bp</p>
                {analysisResults.data.gcContent && (
                  <p className="text-sm">GC Content: {analysisResults.data.gcContent}</p>
                )}
              </div>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">Composition</p>
              <div className="grid grid-cols-5 gap-2">
                {Object.entries(analysisResults.data.compositionPercent).map(([base, percent]) => (
                  <div key={base} className="text-center">
                    <div className="text-sm font-medium">{base}</div>
                    <div className="text-xs">{percent}</div>
                  </div>
                ))}
              </div>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm">
              <p className="font-medium mb-1">Sequences</p>
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b">
                    <th className="text-left pb-1">Name</th>
                    <th className="text-right pb-1">Length</th>
                    <th className="text-left pb-1">Preview</th>
                  </tr>
                </thead>
                <tbody>
                  {analysisResults.data.sequences.map((seq, index) => (
                    <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                      <td className="py-1 truncate" style={{maxWidth: "150px"}}>{seq.header}</td>
                      <td className="py-1 text-right">{seq.length}</td>
                      <td className="py-1 font-mono text-xs truncate">{seq.preview}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
        
        {/* GC Content Analysis */}
        {analysisResults.type === "gc_content" && (
          <div>
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">GC Content Overview</p>
              <p className="text-lg font-bold text-center my-2">{analysisResults.data.overallGcContent}</p>
              <p className="text-sm text-center">Overall GC content across {analysisResults.data.sequenceCount} sequences</p>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm">
              <p className="font-medium mb-1">GC Content by Sequence</p>
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b">
                    <th className="text-left pb-1">Sequence</th>
                    <th className="text-right pb-1">Length</th>
                    <th className="text-right pb-1">GC Content</th>
                  </tr>
                </thead>
                <tbody>
                  {analysisResults.data.sequenceGcContents.map((seq, index) => (
                    <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                      <td className="py-1 truncate" style={{maxWidth: "200px"}}>{seq.header}</td>
                      <td className="py-1 text-right">{seq.length}</td>
                      <td className="py-1 text-right font-medium">{seq.gcContent}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
        
        {/* Protein Translation */}
        {analysisResults.type === "translate" && (
          <div>
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">Translation Results</p>
              <p className="text-sm">Translated {analysisResults.data.sequenceCount} sequences</p>
              
              <div className="mt-3">
                <p className="font-medium text-sm">Translated Sequences:</p>
                <table className="w-full text-sm mt-1">
                  <thead>
                    <tr className="border-b">
                      <th className="text-left pb-1">Sequence</th>
                      <th className="text-right pb-1">Nucleotide Length</th>
                      <th className="text-right pb-1">Protein Length</th>
                    </tr>
                  </thead>
                  <tbody>
                    {analysisResults.data.translatedSequences.map((seq, index) => (
                      <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                        <td className="py-1 truncate" style={{maxWidth: "200px"}}>{seq.header}</td>
                        <td className="py-1 text-right">{seq.nucleotideLength}</td>
                        <td className="py-1 text-right">{seq.proteinLength}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm">
              <p className="font-medium mb-1">FASTA Output</p>
              <div className="font-mono text-xs whitespace-pre-wrap border p-2 rounded max-h-40 overflow-y-auto">
                {analysisResults.data.fastaFormat}
              </div>
              <button
                className="mt-2 text-blue-600 text-sm flex items-center"
                onClick={() => {
                  // Create a download for the FASTA output
                  const blob = new Blob([analysisResults.data.fastaFormat], { type: 'text/plain' });
                  const url = URL.createObjectURL(blob);
                  const a = document.createElement('a');
                  a.href = url;
                  a.download = "translated_sequences.fasta";
                  document.body.appendChild(a);
                  a.click();
                  document.body.removeChild(a);
                  URL.revokeObjectURL(url);
                }}
              >
                <Download className="w-4 h-4 mr-1" />
                Download FASTA
              </button>
            </div>
          </div>
        )}
        
        {/* Feature Statistics */}
        {analysisResults.type === "feature_stats" && (
          <div>
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">Feature Overview</p>
              <p className="text-sm">Total Features: {analysisResults.data.featureCount}</p>
              {analysisResults.data.metadata && (
                <>
                  {analysisResults.data.metadata.name && (
                    <p className="text-sm">Sequence: {analysisResults.data.metadata.name}</p>
                  )}
                  {analysisResults.data.metadata.definition && (
                    <p className="text-sm">Definition: {analysisResults.data.metadata.definition}</p>
                  )}
                </>
              )}
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">Feature Types</p>
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b">
                    <th className="text-left pb-1">Type</th>
                    <th className="text-right pb-1">Count</th>
                  </tr>
                </thead>
                <tbody>
                  {analysisResults.data.featureTypes.map((feature, index) => (
                    <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                      <td className="py-1">{feature.type}</td>
                      <td className="py-1 text-right">{feature.count}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm">
              <p className="font-medium mb-1">Sample Features</p>
              <div className="space-y-2 mt-1">
                {analysisResults.data.topFeatures.map((feature, index) => (
                  <div key={index} className={`p-2 rounded ${index % 2 === 0 ? 'bg-gray-50' : ''}`}>
                    <p className="font-medium text-sm">{feature.type}</p>
                    <p className="text-xs">Location: {feature.location}</p>
                    {Object.keys(feature.qualifiers).length > 0 && (
                      <div className="mt-1">
                        <p className="text-xs font-medium">Qualifiers:</p>
                        <ul className="text-xs list-disc list-inside">
                          {Object.entries(feature.qualifiers).map(([key, value]) => (
                            <li key={key}>{key}: {value}</li>
                          ))}
                        </ul>
                      </div>
                    )}
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}
        
        {/* Structure Analysis */}
        {analysisResults.type === "structure_analysis" && (
          <div>
            <div className="grid grid-cols-2 gap-4 mb-4">
              <div className="bg-white p-3 rounded shadow-sm">
                <p className="font-medium mb-1">Overview</p>
                <p className="text-sm">Title: {analysisResults.data.summary.title}</p>
                <p className="text-sm">Atoms: {analysisResults.data.summary.atomCount}</p>
                <p className="text-sm">HETATM: {analysisResults.data.summary.hetatmCount}</p>
                <p className="text-sm">Chains: {analysisResults.data.summary.chainCount}</p>
              </div>
              
              <div className="bg-white p-3 rounded shadow-sm">
                <p className="font-medium mb-1">Secondary Structure</p>
                <p className="text-sm">Helices: {analysisResults.data.summary.secondaryStructure.helixCount}</p>
                <p className="text-sm">Sheets: {analysisResults.data.summary.secondaryStructure.sheetCount}</p>
                <p className="text-sm">Residue Types: {analysisResults.data.summary.residueTypeCount}</p>
              </div>
            </div>
            
            <div className="grid grid-cols-2 gap-4">
              <div className="bg-white p-3 rounded shadow-sm">
                <p className="font-medium mb-1">Residue Composition</p>
                <div className="max-h-40 overflow-y-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b">
                        <th className="text-left pb-1">Residue</th>
                        <th className="text-right pb-1">Count</th>
                      </tr>
                    </thead>
                    <tbody>
                      {analysisResults.data.residueComposition.map((residue, index) => (
                        <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                          <td className="py-1">{residue.residue}</td>
                          <td className="py-1 text-right">{residue.count}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
              
              <div className="bg-white p-3 rounded shadow-sm">
                <p className="font-medium mb-1">Chain Information</p>
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b">
                      <th className="text-left pb-1">Chain</th>
                      <th className="text-left pb-1">Description</th>
                    </tr>
                  </thead>
                  <tbody>
                    {analysisResults.data.chains.map((chain, index) => (
                      <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                        <td className="py-1">{chain.id}</td>
                        <td className="py-1">{chain.description}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        )}
        
        {/* Conservation Analysis */}
        {analysisResults.type === "conservation" && (
          <div>
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">Conservation Overview</p>
              <p className="text-sm">Sequences: {analysisResults.data.sequenceCount}</p>
              <p className="text-sm">Alignment Length: {analysisResults.data.alignmentLength} positions</p>
              <p className="text-sm">Conservation Score: {analysisResults.data.overallConservation}%</p>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm mb-4">
              <p className="font-medium mb-1">Conserved Regions</p>
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b">
                    <th className="text-left pb-1">Region</th>
                    <th className="text-center pb-1">Length</th>
                    <th className="text-right pb-1">Conservation</th>
                  </tr>
                </thead>
                <tbody>
                  {analysisResults.data.conservedRegions.map((region, index) => (
                    <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : ''}>
                      <td className="py-1">Position {region.start}-{region.end}</td>
                      <td className="py-1 text-center">{region.length}</td>
                      <td className="py-1 text-right">{region.conservation}%</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            
            <div className="bg-white p-3 rounded shadow-sm">
              <p className="font-medium mb-1">Conservation Visualization</p>
              <div className="font-mono text-xs whitespace-pre-wrap overflow-x-auto mt-2">
                {analysisResults.data.visualization.map((line, index) => (
                  <div key={index}>{line}</div>
                ))}
              </div>
            </div>
          </div>
        )}
      </div>
    );
  };
  
  // Render saved jobs
  const renderSavedJobs = () => {
    if (savedJobs.length === 0) {
      return (
        <div className="text-center p-8">
          <Bookmark className="w-12 h-12 mx-auto text-gray-300" />
          <p className="mt-2 text-gray-500">No saved jobs yet</p>
          <p className="text-sm text-gray-400">Convert files and use the save button to store them here</p>
        </div>
      );
    }
    
    return (
      <div className="space-y-3">
        {savedJobs.map(job => (
          <div key={job.id} className="border rounded p-3 hover:bg-gray-50">
            <div className="flex justify-between items-start">
              <div>
                <div className="font-medium">{job.name}</div>
                <div className="text-sm text-gray-500">
                  {job.formatName} · {new Date(job.savedAt).toLocaleString()}
                </div>
              </div>
              <div className="flex space-x-2">
                <button
                  className="p-1 text-blue-600 hover:bg-blue-50 rounded"
                  onClick={() => loadSavedJob(job)}
                  title="Load"
                >
                  <FileText className="w-4 h-4" />
                </button>
                <button
                  className="p-1 text-red-600 hover:bg-red-50 rounded"
                  onClick={() => deleteSavedJob(job.id)}
                  title="Delete"
                >
                  <AlertCircle className="w-4 h-4" />
                </button>
              </div>
            </div>
          </div>
        ))}
      </div>
    );
  };
  
  // Render job history
  const renderJobHistory = () => {
    if (jobHistory.length === 0) {
      return (
        <div className="text-center p-8">
          <Clock className="w-12 h-12 mx-auto text-gray-300" />
          <p className="mt-2 text-gray-500">No history yet</p>
          <p className="text-sm text-gray-400">Your activity will be shown here</p>
        </div>
      );
    }
    
    return (
      <div className="space-y-2">
        {jobHistory.map(job => (
          <div key={job.id} className="border-b pb-2 last:border-b-0">
            <div className="text-sm">
              {job.type === 'upload' && (
                <div className="flex items-center">
                  <Upload className="w-3 h-3 mr-2 text-blue-500" />
                  <span>Uploaded <strong>{job.fileName}</strong> ({FormatDefinitions[job.format]?.name})</span>
                </div>
              )}
              {job.type === 'conversion' && (
                <div className="flex items-center">
                  <RefreshCw className="w-3 h-3 mr-2 text-green-500" />
                  <span>Converted from <strong>{FormatDefinitions[job.sourceFormat]?.name}</strong> to <strong>{FormatDefinitions[job.targetFormat]?.name}</strong></span>
                </div>
              )}
              {job.type === 'analysis' && (
                <div className="flex items-center">
                  <BarChart2 className="w-3 h-3 mr-2 text-purple-500" />
                  <span>Analyzed <strong>{job.fileName}</strong> with {AnalysisTools[job.format]?.find(t => t.id === job.tool)?.name || job.tool}</span>
                </div>
              )}
              {job.type === 'example' && (
                <div className="flex items-center">
                  <FileText className="w-3 h-3 mr-2 text-blue-500" />
                  <span>Loaded example <strong>{job.fileName}</strong></span>
                </div>
              )}
              {job.type === 'url' && (
                <div className="flex items-center">
                  <Globe className="w-3 h-3 mr-2 text-blue-500" />
                  <span>Loaded from URL: <strong>{job.database}</strong></span>
                </div>
              )}
              {job.type === 'search_result' && (
                <div className="flex items-center">
                  <Search className="w-3 h-3 mr-2 text-blue-500" />
                  <span>Loaded search result: <strong>{job.id}</strong> from {job.database}</span>
                </div>
              )}
              <div className="text-xs text-gray-400 ml-5 mt-1">
                {new Date(job.timestamp).toLocaleString()}
              </div>
            </div>
          </div>
        ))}
      </div>
    );
  };
  
  // Render examples tab
  const renderExamples = () => {
    return (
      <div>
        <h2 className="text-lg font-medium mb-3">Example Files</h2>
        <p className="text-gray-600 mb-4">Try out the converter with these pre-loaded example files:</p>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          {Object.entries(ExampleFiles).map(([format, example]) => (
            <div 
              key={format} 
              className={`border rounded p-3 hover:border-blue-300 hover:bg-blue-50 cursor-pointer ${selectedExample === format ? 'border-blue-300 bg-blue-50' : ''}`}
              onClick={() => loadExampleFile(format)}
            >
              <div className="flex items-start">
                {format === 'fasta' && <FileText className="w-6 h-6 mr-2 text-blue-500" />}
                {format === 'genbank' && <FileCode className="w-6 h-6 mr-2 text-green-500" />}
                {format === 'pdb' && <FileCode className="w-6 h-6 mr-2 text-purple-500" />}
                {format === 'clustal' && <List className="w-6 h-6 mr-2 text-orange-500" />}
                <div>
                  <h3 className="font-medium">{example.name}</h3>
                  <p className="text-sm text-gray-600">{example.description}</p>
                  <p className="text-xs text-gray-500 mt-1">{FormatDefinitions[format].name} format</p>
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>
    );
  };
  
  // Render help tab
  const renderHelp = () => {
    return (
      <div>
        <h2 className="text-lg font-medium mb-3">Help & Information</h2>
        
        <div className="mb-4 p-4 bg-blue-50 rounded border border-blue-200">
          <h3 className="font-medium text-blue-800 mb-2 flex items-center">
            <HelpCircle className="w-5 h-5 mr-2" />
            About This Tool
          </h3>
          <p className="text-sm text-blue-800 mb-2">
            The Biology Format Converter allows you to convert between different file formats commonly used in bioinformatics and computational biology. 
            All processing happens in your browser - no data is sent to any server.
          </p>
          <p className="text-sm text-blue-800">
            This tool is designed for researchers, students, and professionals working with biological data who need to quickly convert between different file formats.
          </p>
        </div>
        
        <div className="space-y-4">
          <div className="border-b pb-2">
            <h3 className="font-medium mb-2">How to Use</h3>
            <ol className="list-decimal list-inside text-sm space-y-2">
              <li>Upload a file using the file dialog or drag-and-drop</li>
              <li>The tool will automatically detect the format of your file</li>
              <li>Select the target format you want to convert to</li>
              <li>Adjust any conversion options if needed</li>
              <li>Click "Convert" to process your file</li>
              <li>Preview the result and download the converted file</li>
            </ol>
          </div>
          
          <div className="border-b pb-2">
            <h3 className="font-medium mb-2">Supported Formats</h3>
            <div className="text-sm space-y-3">
              {Object.entries(FormatCategories).map(([category, info]) => (
                <div key={category}>
                  <p className="font-medium">{info.name}</p>
                  <ul className="list-disc list-inside text-sm">
                    {Object.entries(FormatDefinitions)
                      .filter(([_, formatInfo]) => formatInfo.category === category)
                      .map(([key, formatInfo]) => (
                        <li key={key}>{formatInfo.name} ({formatInfo.extensions.join(', ')})</li>
                      ))
                    }
                  </ul>
                </div>
              ))}
            </div>
          </div>
          
          <div>
            <h3 className="font-medium mb-2">Analysis Tools</h3>
            <ul className="list-disc list-inside text-sm space-y-1 text-gray-700">
              <li>Sequence Statistics: Get basic statistics about sequences</li>
              <li>GC Content Analysis: Calculate GC content and distribution</li>
              <li>Translation: Convert nucleotide sequences to protein</li>
              <li>Feature Analysis: Analyze annotated features in GenBank files</li>
              <li>Structure Analysis: Analyze protein structures from PDB files</li>
              <li>Conservation Analysis: Analyze multiple sequence alignments</li>
            </ul>
          </div>
        </div>
      </div>
    );
  };
  
  // Main render
  return (
    <div className="max-w-4xl mx-auto p-4">
      <div className="bg-white rounded-lg shadow p-6">
        <h1 className="text-2xl font-bold mb-1">Biology Format Converter</h1>
        <p className="text-gray-600 mb-4">Convert between common biological data formats</p>
        
        <div className="flex border-b mb-4">
          <button
            className={`px-4 py-2 font-medium ${activeTab === 'upload' ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
            onClick={() => setActiveTab('upload')}
          >
            Upload File
          </button>
          <button
            className={`px-4 py-2 font-medium ${activeTab === 'url' ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
            onClick={() => setActiveTab('url')}
          >
            URL / Database
          </button>
          <button
            className={`px-4 py-2 font-medium ${activeTab === 'examples' ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
            onClick={() => setActiveTab('examples')}
          >
            Examples
          </button>
          <button
            className={`px-4 py-2 font-medium ${activeTab === 'saved' ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
            onClick={() => setActiveTab('saved')}
          >
            Saved Jobs
          </button>
          <button
            className={`px-4 py-2 font-medium ${activeTab === 'history' ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
            onClick={() => setActiveTab('history')}
          >
            History
          </button>
          <button
            className={`px-4 py-2 font-medium ${activeTab === 'help' ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
            onClick={() => setActiveTab('help')}
          >
            Help
          </button>
        </div>
        
        {/* Upload Tab */}
        {activeTab === 'upload' && (
          <>
            <div 
              className={`border-2 ${dragActive ? 'border-blue-400 bg-blue-50' : 'border-dashed border-gray-300 bg-gray-50'} rounded-lg p-8 text-center`}
              onDragEnter={handleDrag}
              onDragOver={handleDrag}
              onDragLeave={handleDrag}
              onDrop={handleDrop}
            >
              <input
                type="file"
                ref={fileInputRef}
                onChange={handleFileUpload}
                className="hidden"
              />
              
              {!uploadedFile ? (
                <div>
                  <Upload className="w-12 h-12 mx-auto text-gray-400" />
                  <p className="mt-2 text-gray-600">Drag and drop your file here or</p>
                  <button
                    className="mt-2 bg-blue-600 text-white px-4 py-2 rounded"
                    onClick={() => fileInputRef.current.click()}
                  >
                    Browse Files
                  </button>
                  
                  <div className="mt-4">
                    <button
                      className="text-blue-600 text-sm flex items-center"
                      onClick={() => setShowFormatInfo(!showFormatInfo)}
                    >
                      <ChevronDown className={`w-4 h-4 mr-1 transition-transform ${showFormatInfo ? 'rotate-180' : ''}`} />
                      {showFormatInfo ? 'Hide supported formats' : 'Show supported formats'}
                    </button>
                    
                    {showFormatInfo && (
                      <div className="mt-2 text-sm border rounded p-3 bg-gray-50">
                        <h3 className="font-medium mb-2">Supported Formats:</h3>
                        
                        {Object.entries(FormatCategories).map(([category, info]) => (
                          <div key={category} className="mb-3">
                            <h4 className="font-medium text-gray-700">{info.name}</h4>
                            <p className="text-xs text-gray-600 mb-1">{info.description}</p>
                            
                            <div className="space-y-2">
                              {Object.entries(FormatDefinitions)
                                .filter(([_, formatInfo]) => formatInfo.category === category)
                                .map(([key, format]) => (
                                  <div key={key} className="border-b pb-2">
                                    <p className="font-medium">{format.name} ({format.extensions.join(', ')})</p>
                                    <p className="text-gray-600">{format.description}</p>
                                  </div>
                                ))}
                            </div>
                          </div>
                        ))}
                        
                        <h3 className="font-medium mt-4 mb-2">Supported Conversions:</h3>
                        <ul className="list-disc list-inside">
                          {Object.entries(ConversionMatrix).flatMap(([from, toFormats]) => 
                            toFormats.map(to => (
                              <li key={`${from}-${to}`}>
                                {FormatDefinitions[from].name} → {FormatDefinitions[to].name}
                              </li>
                            ))
                          )}
                        </ul>
                      </div>
                    )}
                  </div>
                </div>
              ) : (
                <div>
                  <div className="flex items-center justify-center mb-4">
                    <FileText className="w-8 h-8 text-blue-600 mr-2" />
                    <div className="text-left">
                      <p className="font-medium">{uploadedFile.name}</p>
                      <p className="text-sm text-gray-500">
                        {(uploadedFile.size / 1024).toFixed(1)} KB · 
                        {detectedFormat ? 
                          ` Detected as ${FormatDefinitions[detectedFormat].name}` : 
                          ' Format not detected'}
                      </p>
                    </div>
                  </div>
                  
                  <div className="flex space-x-2 justify-center">
                    <button
                      className="text-sm text-blue-600"
                      onClick={resetState}
                    >
                      Upload a different file
                    </button>
                    
                    {fileType && (
                      <button
                        className="text-sm text-blue-600 flex items-center"
                        onClick={saveCurrentJob}
                      >
                        <Save className="w-4 h-4 mr-1" />
                        Save
                      </button>
                    )}
                  </div>
                </div>
              )}
            </div>
            
            {error && (
              <div className="mt-4 bg-red-100 border border-red-200 text-red-700 p-3 rounded flex items-start">
                <AlertCircle className="w-5 h-5 mr-2 flex-shrink-0 mt-0.5" />
                <p>{error}</p>
              </div>
            )}
            
            {processingProgress > 0 && processingProgress < 100 && (
              <div className="mt-4">
                <div className="w-full bg-gray-200 rounded-full h-2.5">
                  <div className="bg-blue-600 h-2.5 rounded-full" style={{ width: `${processingProgress}%` }}></div>
                </div>
                <p className="text-sm text-gray-600 mt-1">Processing file... {processingProgress}%</p>
              </div>
            )}
            
            {uploadedFile && (
              <div>
                {renderFormatOptions()}
                {renderAnalysisTools()}
                {renderPreview()}
                {renderConversionResult()}
                {renderAnalysisResults()}
              </div>
            )}
          </>
        )}
        
        {/* URL Tab */}
        {activeTab === 'url' && (
          <UrlLoader 
            urlInput={urlInput}
            setUrlInput={setUrlInput}
            isLoadingUrl={isLoadingUrl}
            processUrlInput={processUrlInput}
            databaseSearch={databaseSearch}
            setDatabaseSearch={setDatabaseSearch}
            isSearching={isSearching}
            performDatabaseSearch={performDatabaseSearch}
            searchResults={searchResults}
            loadSearchResult={loadSearchResult}
            databaseResources={DatabaseResources}
            error={error}
            setError={setError}
          />
        )}
        
        {/* Examples Tab */}
        {activeTab === 'examples' && renderExamples()}
        
        {/* Saved Jobs Tab */}
        {activeTab === 'saved' && (
          <div>
            <h2 className="text-lg font-medium mb-3">Saved Jobs</h2>
            {renderSavedJobs()}
          </div>
        )}
        
        {/* History Tab */}
        {activeTab === 'history' && (
          <div>
            <h2 className="text-lg font-medium mb-3">Activity History</h2>
            {renderJobHistory()}
          </div>
        )}
        
        {/* Help Tab */}
        {activeTab === 'help' && renderHelp()}
        
        <div className="mt-6 text-xs text-gray-500 border-t pt-4">
          <p>All processing is done in your browser - no data is sent to a server.</p>
          <p>Created for scientific research. Designed to handle commonly used biological data formats.</p>
        </div>
      </div>
      
      {/* Notification */}
      {showNotification && (
        <div className="fixed bottom-4 right-4 bg-blue-600 text-white px-4 py-2 rounded-md shadow-lg">
          {notificationMessage}
        </div>
      )}
    </div>
  );
};

// UIComponents.js
import React from 'react';

export const FileDropZone = ({ children, onDrop, onDragOver, onDragLeave, dragActive }) => {
  return (
    <div 
      className={`border-2 ${dragActive ? 'border-blue-400 bg-blue-50' : 'border-dashed border-gray-300 bg-gray-50'} rounded-lg p-8 text-center`}
      onDragOver={onDragOver}
      onDragLeave={onDragLeave}
      onDrop={onDrop}
    >
      {children}
    </div>
  );
};

export const Notification = ({ message, onClose }) => {
  return (
    <div className="fixed bottom-4 right-4 bg-blue-600 text-white px-4 py-2 rounded-md shadow-lg flex justify-between items-center">
      <span>{message}</span>
      {onClose && (
        <button className="ml-3 text-white" onClick={onClose}>×</button>
      )}
    </div>
  );
};

export const TabNavigation = ({ tabs, activeTab, onTabChange }) => {
  return (
    <div className="flex border-b mb-4">
      {tabs.map(tab => (
        <button
          key={tab.id}
          className={`px-4 py-2 font-medium ${activeTab === tab.id ? 'border-b-2 border-blue-500 text-blue-600' : 'text-gray-600'}`}
          onClick={() => onTabChange(tab.id)}
        >
          {tab.label}
        </button>
      ))}
    </div>
  );
};

export const ConversionOptions = ({ options, setOptions, sourceFormat, targetFormat }) => {
  return (
    <div className="mt-3 border rounded-md p-3 bg-gray-50">
      <h3 className="text-sm font-medium mb-2">Conversion Options</h3>
      
      <div className="grid grid-cols-2 gap-2">
        {(sourceFormat === 'fasta' || targetFormat === 'fasta') && (
          <label className="flex items-center">
            <input
              type="checkbox"
              checked={options.formatSequence}
              onChange={(e) => setOptions({...options, formatSequence: e.target.checked})}
              className="mr-2"
            />
            <span className="text-sm">Format sequence with line breaks</span>
          </label>
        )}
        
        {sourceFormat === 'fastq' && (
          <label className="flex items-center">
            <input
              type="checkbox"
              checked={options.includeQuality}
              onChange={(e) => setOptions({...options, includeQuality: e.target.checked})}
              className="mr-2"
            />
            <span className="text-sm">Include quality scores</span>
          </label>
        )}
        
        {(sourceFormat === 'clustal' || sourceFormat === 'stockholm') && (
          <label className="flex items-center">
            <input
              type="checkbox"
              checked={options.removeGaps}
              onChange={(e) => setOptions({...options, removeGaps: e.target.checked})}
              className="mr-2"
            />
            <span className="text-sm">Remove alignment gaps</span>
          </label>
        )}
        
        {/* Add more format-specific options here */}
      </div>
    </div>
  );
};


// Process URL input
const processUrlInput = (url) => {
  if (!url.trim()) {
    setError("Please enter a URL");
    return;
  }
  
  setIsLoadingUrl(true);
  setError(null);
  
  // Detect if it's a database URL
  const dbInfo = getDatabaseFromUrl(url);
  
  fetch(url)
    .then(response => {
      if (!response.ok) {
        throw new Error(`HTTP error! Status: ${response.status}`);
      }
      return response.text();
    })
    .then(content => {
      // Auto-detect format
      const detectedFormat = dbInfo.database ? 
        DatabaseUrlToFormat[dbInfo.database] || detectFormatByContent(content) : 
        detectFormatByContent(content);
      
      if (!detectedFormat) {
        throw new Error("Could not detect file format");
      }
      
      // Create file object
      const fileName = url.split('/').pop() || 'downloaded_file';
      const file = new File([content], fileName, { type: 'text/plain' });
      
      setUploadedFile(file);
      setFileContent(content);
      setFileType(detectedFormat);
      setDetectedFormat(detectedFormat);
      generatePreview(content, detectedFormat);
      
      // Add to job history
      addToHistory({
        type: 'url',
        format: detectedFormat,
        url: url,
        database: dbInfo.database || 'unknown',
        timestamp: new Date().toISOString()
      });
      
      showNotificationMessage(`Loaded file from URL: ${fileName}`);
      setIsLoadingUrl(false);
      setActiveTab('upload'); // Switch to upload tab to show preview
    })
    .catch(err => {
      setError(`Error loading URL: ${err.message}`);
      setIsLoadingUrl(false);
    });
};

// Search database
const performDatabaseSearch = (database, subDatabase, searchTerm) => {
  if (!searchTerm.trim()) {
    setError("Please enter a search term");
    return;
  }
  
  setIsSearching(true);
  setError(null);
  setSearchResults([]);
  
  // Generate search URL
  const searchUrl = generateSearchUrl(database, subDatabase, searchTerm);
  
  if (!searchUrl) {
    setError("Could not generate search URL");
    setIsSearching(false);
    return;
  }
  
  // Perform search
  searchDatabase(
    database,
    subDatabase,
    searchTerm,
    (results) => {
      setSearchResults(results.results || []);
      setIsSearching(false);
      
      if (results.results.length === 0) {
        setError(`No results found for "${searchTerm}"`);
      }
    },
    (error) => {
      setError(`Search error: ${error}`);
      setIsSearching(false);
    }
  );
};

// Load search result
const loadSearchResult = (result) => {
  setIsLoadingUrl(true);
  setError(null);
  
  const database = result.database;
  const subDatabase = result.subDatabase;
  const accessionId = result.id;
  
  // Generate URL for fetching the entry
  const fetchUrl = generateDatabaseUrl(database, subDatabase, accessionId);
  
  if (!fetchUrl) {
    setError("Could not generate URL for the selected entry");
    setIsLoadingUrl(false);
    return;
  }
  
  // Fetch the entry
  fetchFromUrl(
    fetchUrl,
    (result) => {
      const content = result.content;
      const format = result.format;
      
      if (!content || !format) {
        setError("Could not retrieve data from the selected entry");
        setIsLoadingUrl(false);
        return;
      }
      
      // Create file object
      const fileName = `${accessionId}.${format}`;
      const file = new File([content], fileName, { type: 'text/plain' });
      
      setUploadedFile(file);
      setFileContent(content);
      setFileType(format);
      setDetectedFormat(format);
      generatePreview(content, format);
      
      // Add to job history
      addToHistory({
        type: 'search_result',
        format: format,
        id: accessionId,
        database: database,
        timestamp: new Date().toISOString()
      });
      
      showNotificationMessage(`Loaded ${accessionId} from ${database}`);
      setIsLoadingUrl(false);
      setActiveTab('upload'); // Switch to upload tab to show preview
    },
    (error) => {
      setError(`Error loading entry: ${error}`);
      setIsLoadingUrl(false);
    }
  );
};


  return (
    <div className="max-w-4xl mx-auto p-4">
      <div className="bg-white rounded-lg shadow p-6">
        <h1 className="text-2xl font-bold mb-1">Biology Format Converter</h1>
        <p className="text-gray-600 mb-4">Convert between common biological data formats</p>
        
        {/* Error display with details for advanced users */}
        {error && (
          <div className="mt-4 bg-red-100 border border-red-200 text-red-700 p-3 rounded flex items-start">
            <AlertCircle className="w-5 h-5 mr-2 flex-shrink-0 mt-0.5" />
            <div>
              <p className="font-medium">{error}</p>
              {errorDetails && (
                <details className="mt-1 text-xs">
                  <summary className="cursor-pointer">Technical details</summary>
                  <p className="mt-1">Source: {errorDetails.source || 'Unknown'}</p>
                  {errorDetails.message && <p>Message: {errorDetails.message}</p>}
                  {errorDetails.stack && (
                    <pre className="mt-1 whitespace-pre-wrap overflow-x-auto max-h-40">
                      {errorDetails.stack}
                    </pre>
                  )}
                </details>
              )}
            </div>
          </div>
        )}
        
        {/* Rest of the UI components */}
        {/* ... */}
        
        {/* Notification */}
        {showNotification && (
          <div className="fixed bottom-4 right-4 bg-blue-600 text-white px-4 py-2 rounded-md shadow-lg">
            {notificationMessage}
          </div>
        )}
      </div>
    </div>
  );
};

export default BioFormatConverter;
