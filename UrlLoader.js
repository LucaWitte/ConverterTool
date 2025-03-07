// UrlLoader.js - Component for loading data from URLs and database searches

import React, { useState } from 'react';
import { Globe, Search, Database, Download, Link, ChevronDown, List, AlertCircle, RefreshCw } from 'lucide-react';

/**
 * URL Loader component
 * @param {object} props - Component props
 * @returns {JSX.Element} React component
 */
export const UrlLoader = ({
  urlInput,
  setUrlInput,
  isLoadingUrl,
  processUrlInput,
  databaseSearch,
  setDatabaseSearch,
  isSearching,
  performDatabaseSearch,
  searchResults,
  loadSearchResult,
  databaseResources,
  error,
  setError
}) => {
  const [selectedDatabase, setSelectedDatabase] = useState('ncbi');
  const [selectedSubDatabase, setSelectedSubDatabase] = useState('nucleotide');
  const [showDbSelector, setShowDbSelector] = useState(false);
  
  // Handle URL input change
  const handleUrlChange = (e) => {
    setUrlInput(e.target.value);
    if (error) setError(null);
  };
  
  // Handle database selection
  const handleDatabaseSelect = (dbId) => {
    setSelectedDatabase(dbId);
    // Select first sub-database by default
    const firstSubDb = databaseResources[dbId].databases[0].id;
    setSelectedSubDatabase(firstSubDb);
  };
  
  // Handle sub-database selection
  const handleSubDatabaseSelect = (subDbId) => {
    setSelectedSubDatabase(subDbId);
  };
  
  // Handle search input change
  const handleSearchChange = (e) => {
    setDatabaseSearch(e.target.value);
    if (error) setError(null);
  };
  
  // Handle search submission
  const handleSearchSubmit = (e) => {
    e.preventDefault();
    if (databaseSearch.trim()) {
      performDatabaseSearch(selectedDatabase, selectedSubDatabase, databaseSearch.trim());
    }
  };
  
  // Get current database and sub-database
  const currentDatabase = databaseResources[selectedDatabase];
  const currentSubDatabase = currentDatabase?.databases.find(db => db.id === selectedSubDatabase);
  
  return (
    <div>
      <h2 className="text-lg font-medium mb-3">Load from URL or Database</h2>
      
      {/* URL input section */}
      <div className="mb-5">
        <div className="text-gray-600 mb-2 flex items-center">
          <Globe className="w-4 h-4 mr-2" />
          <h3>Direct URL</h3>
        </div>
        
        <div className="flex">
          <input
            type="text"
            className="flex-1 p-2 border rounded-l focus:outline-none focus:border-blue-500"
            placeholder="Enter URL to a biological data file..."
            value={urlInput}
            onChange={handleUrlChange}
          />
          <button
            className="bg-blue-600 text-white px-4 py-2 rounded-r flex items-center"
            onClick={() => processUrlInput(urlInput)}
            disabled={isLoadingUrl || !urlInput.trim()}
          >
            {isLoadingUrl ? (
              <>
                <RefreshCw className="w-4 h-4 mr-2 animate-spin" />
                Loading...
              </>
            ) : (
              <>
                <Download className="w-4 h-4 mr-2" />
                Load
              </>
            )}
          </button>
        </div>
        
        <div className="mt-1 text-sm text-gray-500">
          Enter a direct URL to a file in any supported format, or a database entry URL.
        </div>
      </div>
      
      {/* Database search section */}
      <div className="mb-5">
        <div className="text-gray-600 mb-2 flex items-center">
          <Database className="w-4 h-4 mr-2" />
          <h3>Database Search</h3>
        </div>
        
        <div className="rounded border p-3 bg-gray-50">
          <div className="mb-3">
            <div className="text-sm font-medium mb-1">Select Database:</div>
            <div className="relative">
              <button
                className="w-full p-2 border rounded bg-white flex justify-between items-center"
                onClick={() => setShowDbSelector(!showDbSelector)}
              >
                <span>{currentDatabase?.name || 'Select Database'}</span>
                <ChevronDown className="w-4 h-4" />
              </button>
              
              {showDbSelector && (
                <div className="absolute left-0 right-0 mt-1 bg-white border rounded shadow-lg z-10 max-h-64 overflow-y-auto">
                  {Object.entries(databaseResources).map(([dbId, db]) => (
                    <button
                      key={dbId}
                      className={`w-full text-left p-2 hover:bg-blue-50 ${selectedDatabase === dbId ? 'bg-blue-100' : ''}`}
                      onClick={() => {
                        handleDatabaseSelect(dbId);
                        setShowDbSelector(false);
                      }}
                    >
                      <div className="font-medium">{db.name}</div>
                      <div className="text-xs text-gray-600">{db.description}</div>
                    </button>
                  ))}
                </div>
              )}
            </div>
          </div>
          
          <div className="mb-3">
            <div className="text-sm font-medium mb-1">Select Category:</div>
            <div className="grid grid-cols-2 gap-2">
              {currentDatabase?.databases.map(subDb => (
                <button
                  key={subDb.id}
                  className={`p-2 border rounded text-left ${selectedSubDatabase === subDb.id ? 'bg-blue-100 border-blue-300' : 'bg-white'}`}
                  onClick={() => handleSubDatabaseSelect(subDb.id)}
                >
                  <div className="font-medium">{subDb.name}</div>
                  <div className="text-xs text-gray-600">{subDb.description}</div>
                </button>
              ))}
            </div>
          </div>
          
          <form onSubmit={handleSearchSubmit}>
            <div className="flex">
              <input
                type="text"
                className="flex-1 p-2 border rounded-l focus:outline-none focus:border-blue-500"
                placeholder={`Search ${currentDatabase?.name} ${currentSubDatabase?.name}...`}
                value={databaseSearch}
                onChange={handleSearchChange}
              />
              <button
                type="submit"
                className="bg-blue-600 text-white px-4 py-2 rounded-r flex items-center"
                disabled={isSearching || !databaseSearch.trim()}
              >
                {isSearching ? (
                  <>
                    <RefreshCw className="w-4 h-4 mr-2 animate-spin" />
                    Searching...
                  </>
                ) : (
                  <>
                    <Search className="w-4 h-4 mr-2" />
                    Search
                  </>
                )}
              </button>
            </div>
          </form>
          
          <div className="mt-1 text-xs text-gray-500">
            {currentSubDatabase?.exampleIds?.length > 0 && (
              <div>
                Examples: {currentSubDatabase.exampleIds.map((id, index) => (
                  <button
                    key={id}
                    className="text-blue-600 hover:underline"
                    onClick={() => setDatabaseSearch(id)}
                  >
                    {id}{index < currentSubDatabase.exampleIds.length - 1 ? ', ' : ''}
                  </button>
                ))}
              </div>
            )}
          </div>
        </div>
      </div>
      
      {/* Error message */}
      {error && (
        <div className="bg-red-100 border border-red-200 text-red-700 p-3 rounded flex items-start my-3">
          <AlertCircle className="w-5 h-5 mr-2 flex-shrink-0 mt-0.5" />
          <p>{error}</p>
        </div>
      )}
      
      {/* Search results */}
      {searchResults.length > 0 && (
        <div className="mt-4">
          <h3 className="text-lg font-medium mb-2">Search Results</h3>
          <div className="border rounded divide-y">
            {searchResults.map((result, index) => (
              <div key={index} className="p-3 hover:bg-gray-50">
                <div className="flex justify-between items-start">
                  <div>
                    <div className="font-medium">{result.title || result.id}</div>
                    {result.description && (
                      <div className="text-sm text-gray-600">{result.description}</div>
                    )}
                    {result.organism && (
                      <div className="text-xs text-gray-500">{result.organism}</div>
                    )}
                    <div className="text-xs text-gray-500 mt-1">
                      ID: {result.id} · Database: {result.database} · Type: {result.subDatabase}
                    </div>
                  </div>
                  <button
                    className="bg-blue-600 text-white px-3 py-1 rounded text-sm flex items-center"
                    onClick={() => loadSearchResult(result)}
                  >
                    <Link className="w-3 h-3 mr-1" />
                    Load
                  </button>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
};

export default UrlLoader;
