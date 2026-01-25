import React, { useState, useEffect } from 'react';
import { Card, Button, Form, Row, Col, Alert, Badge, Spinner } from 'react-bootstrap';
import ReactPaginate from 'react-paginate';
import Select from 'react-select';
import MoleculeViewer from './MoleculeViewer';

function MoleculeReviewer({ onReviewUpdate }) {
  const [molecules, setMolecules] = useState([]);
  const [loading, setLoading] = useState(true);
  const [currentPage, setCurrentPage] = useState(1);
  const [totalPages, setTotalPages] = useState(0);
  const [filters, setFilters] = useState({
    dataset: 'all',
    reviewStatus: 'all',
    search: '',
    prioritizeMisclassified: true
  });
  const [reviewData, setReviewData] = useState({});
  const [submitStatus, setSubmitStatus] = useState({});

  const datasetOptions = [
    { value: 'all', label: 'All Datasets' },
    { value: 'oecd', label: 'OECD' },
    { value: 'enhanced', label: 'Enhanced' },
    { value: 'timing', label: 'Timing' },
    { value: 'complex_branched', label: 'Complex Branched' },
    { value: 'non_fluorinated', label: 'Non-Fluorinated' }
  ];

  const reviewStatusOptions = [
    { value: 'all', label: 'All Molecules' },
    { value: 'unreviewed', label: 'Unreviewed Only' },
    { value: 'reviewed', label: 'Reviewed Only' }
  ];

  useEffect(() => {
    fetchMolecules();
  }, [currentPage, filters]);

  const fetchMolecules = async () => {
    setLoading(true);
    try {
      const params = new URLSearchParams({
        page: currentPage,
        limit: 10,
        dataset: filters.dataset,
        reviewStatus: filters.reviewStatus,
        search: filters.search,
        prioritizeMisclassified: filters.prioritizeMisclassified
      });

      const response = await fetch(`/api/molecules?${params}`);
      const data = await response.json();
      
      setMolecules(data.molecules);
      setTotalPages(data.pagination.totalPages);
      
      // Initialize review data for new molecules
      const newReviewData = {};
      data.molecules.forEach(mol => {
        newReviewData[mol.id] = {
          pfasgroupsCorrect: mol.pfasgroups_correct,
          atlasCorrect: mol.atlas_correct,
          reviewerNotes: mol.reviewer_notes || '',
          isPfas: mol.manual_is_pfas,
          correctGroups: mol.manual_correct_groups || [],
          correctClassification: mol.manual_correct_classification || ''
        };
      });
      setReviewData(prev => ({ ...prev, ...newReviewData }));
      
    } catch (error) {
      console.error('Error fetching molecules:', error);
    }
    setLoading(false);
  };

  const handleFilterChange = (field, value) => {
    setFilters(prev => ({ ...prev, [field]: value }));
    setCurrentPage(1);
  };

  const handleReviewChange = (moleculeId, field, value) => {
    setReviewData(prev => ({
      ...prev,
      [moleculeId]: {
        ...prev[moleculeId],
        [field]: value
      }
    }));
  };

  const submitReview = async (moleculeId) => {
    setSubmitStatus(prev => ({ ...prev, [moleculeId]: 'submitting' }));
    
    try {
      const reviewDataForMolecule = reviewData[moleculeId];
      
      const response = await fetch('/api/review', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          moleculeId,
          pfasgroupsCorrect: reviewDataForMolecule.pfasgroupsCorrect,
          atlasCorrect: reviewDataForMolecule.atlasCorrect,
          reviewerNotes: reviewDataForMolecule.reviewerNotes,
          reviewerName: 'Manual Reviewer', // Could be made configurable
          isPfas: reviewDataForMolecule.isPfas,
          correctGroups: reviewDataForMolecule.correctGroups,
          correctClassification: reviewDataForMolecule.correctClassification
        }),
      });

      if (response.ok) {
        setSubmitStatus(prev => ({ ...prev, [moleculeId]: 'success' }));
        onReviewUpdate && onReviewUpdate();
        
        // Clear success message after 2 seconds
        setTimeout(() => {
          setSubmitStatus(prev => ({ ...prev, [moleculeId]: null }));
        }, 2000);
      } else {
        setSubmitStatus(prev => ({ ...prev, [moleculeId]: 'error' }));
      }
    } catch (error) {
      console.error('Error submitting review:', error);
      setSubmitStatus(prev => ({ ...prev, [moleculeId]: 'error' }));
    }
  };

  const handlePageChange = (selectedPage) => {
    setCurrentPage(selectedPage.selected + 1);
  };

  if (loading) {
    return <div className="text-center"><Spinner animation="border" /> Loading molecules...</div>;
  }

  return (
    <div>
      <h1>🔬 Molecule Reviewer</h1>

      {/* Legend for component types color coding */}
      <Alert variant="info" className="mb-3">
        <strong>Component Types Legend:</strong> &nbsp;
        <span className="badge bg-primary me-2">per</span> Perfluoroalkyl &nbsp;
        <span className="badge bg-success me-2">poly</span> Polyfluoroalkyl &nbsp;
        <span className="badge bg-info me-2">per,poly</span> Mixed &nbsp;
        <span className="badge bg-warning text-dark me-2">cyc</span> Cyclic &nbsp;
        <span className="badge bg-secondary me-2">—</span> No data
      </Alert>

      {/* Filters */}
      <Card className="filter-section mb-4">
        <Card.Body>
          <h5>Filters</h5>
          <Row>
            <Col md={3}>
              <Form.Group>
                <Form.Label>Dataset</Form.Label>
                <Select
                  value={datasetOptions.find(opt => opt.value === filters.dataset)}
                  onChange={(option) => handleFilterChange('dataset', option.value)}
                  options={datasetOptions}
                />
              </Form.Group>
            </Col>
            <Col md={3}>
              <Form.Group>
                <Form.Label>Review Status</Form.Label>
                <Select
                  value={reviewStatusOptions.find(opt => opt.value === filters.reviewStatus)}
                  onChange={(option) => handleFilterChange('reviewStatus', option.value)}
                  options={reviewStatusOptions}
                />
              </Form.Group>
            </Col>
            <Col md={3}>
              <Form.Group>
                <Form.Label>Search</Form.Label>
                <Form.Control
                  type="text"
                  placeholder="Search by SMILES or group name..."
                  value={filters.search}
                  onChange={(e) => handleFilterChange('search', e.target.value)}
                />
              </Form.Group>
            </Col>
            <Col md={3}>
              <Form.Group>
                <Form.Label>Prioritization</Form.Label>
                <Form.Check
                  type="switch"
                  id="prioritize-switch"
                  label="Prioritize Misclassified"
                  checked={filters.prioritizeMisclassified}
                  onChange={(e) => handleFilterChange('prioritizeMisclassified', e.target.checked)}
                />
                <Form.Text className="text-muted">
                  Show unreviewed misclassified entries first
                </Form.Text>
              </Form.Group>
            </Col>
          </Row>
        </Card.Body>
      </Card>

      {/* Molecules */}
      {molecules.map(molecule => (
        <Card key={molecule.id} className={`molecule-card mb-4 priority-${molecule.priority_category?.replace('_', '-') || 'default'}`}>
          <Card.Body>
            <Row>
              <Col md={6}>
                <div className="molecule-info">
                  <div className="d-flex align-items-center mb-2">
                    <h5 className="mb-0">Molecule #{molecule.id}</h5>
                    {molecule.priority_category === 'misclassified_unreviewed' && (
                      <Badge bg="danger" className="ms-2">⚠️ Priority Review Needed</Badge>
                    )}
                    {molecule.priority_category === 'unreviewed' && (
                      <Badge bg="warning" className="ms-2">📋 Unreviewed</Badge>
                    )}
                    {molecule.priority_category === 'misclassified_reviewed' && (
                      <Badge bg="info" className="ms-2">🔄 Previously Misclassified</Badge>
                    )}
                    {molecule.priority_category === 'reviewed_correct' && (
                      <Badge bg="success" className="ms-2">✅ Reviewed Correct</Badge>
                    )}
                  </div>
                  <p><strong>SMILES:</strong> <code>{molecule.smiles}</code></p>
                  <p><strong>Dataset:</strong> <Badge bg="info">{molecule.dataset_type.toUpperCase()}</Badge></p>
                  {molecule.group_name && <p><strong>Group:</strong> {molecule.group_name}</p>}
                  {molecule.molecular_weight && <p><strong>MW:</strong> {molecule.molecular_weight.toFixed(2)}</p>}
                  {molecule.chain_length && <p><strong>Chain Length:</strong> {molecule.chain_length}</p>}
                </div>

                {/* Molecule Structure */}
                <div className="molecule-structure">
                  <MoleculeViewer smiles={molecule.smiles} />
                </div>
              </Col>

              <Col md={6}>
                {/* Classification Results */}
                <div className="classification-results">
                  {/* PFASGroups Results */}
                  <div className={`algorithm-result ${molecule.pfasgroups_success ? 'success' : 'error'}`}>
                    <h6>🔬 PFASGroups</h6>
                    {molecule.pfasgroups_success ? (
                      <div>
                        <p><strong>Detected Groups:</strong></p>
                        <div>
                          {molecule.pfasgroups_detected && molecule.pfasgroups_detected.map(group => {
                            // Handle both enriched objects {id, name, matchedPathTypeFull} and simple strings/numbers
                            const groupKey = typeof group === 'object' ? group.id : group;
                            const groupName = typeof group === 'object' ? group.name : group;
                            const pathType = typeof group === 'object' ? group.matchedPathTypeFull : null;
                            const pathAbbrev = typeof group === 'object' ? group.matchedPathType : null;
                            
                            // Color code by path type - handle comma-separated types
                            let badgeColor = 'secondary';
                            if (pathType) {
                              if (pathType.includes('Perfluoroalkyl') && pathType.includes('Polyfluoroalkyl')) {
                                badgeColor = 'info'; // Cyan for mixed
                              } else if (pathType.includes('Perfluoroalkyl') || pathType.includes('Perfluoro')) {
                                badgeColor = 'primary'; // Blue for Perfluoroalkyl
                              } else if (pathType.includes('Polyfluoroalkyl') || pathType.includes('Polyfluoro')) {
                                badgeColor = 'success'; // Green for Polyfluoroalkyl
                              } else if (pathType.includes('cyclic')) {
                                badgeColor = 'warning'; // Yellow for cyclic
                              }
                            }
                            
                            return (
                              <Badge 
                                key={groupKey} 
                                bg={badgeColor} 
                                className="me-1 mb-1"
                                title={pathType ? `Components: ${pathType}` : 'No component type data'}
                              >
                                {groupName} {pathAbbrev && <small>({pathAbbrev})</small>}
                              </Badge>
                            );
                          })}
                        </div>
                        {molecule.pfasgroups_detected_definitions && molecule.pfasgroups_detected_definitions.length > 0 && (
                          <>
                            <p className="mt-2"><strong>Detected PFAS Definitions:</strong></p>
                            <div>
                              {molecule.pfasgroups_detected_definitions.map((defId, idx) => (
                                <Badge key={idx} bg="info" className="me-1 mb-1">D{defId}</Badge>
                              ))}
                            </div>
                          </>
                        )}
                        <p><small>Time: {(molecule.pfasgroups_time * 1000).toFixed(2)}ms</small></p>
                      </div>
                    ) : (
                      <p className="text-danger">Error: {molecule.pfasgroups_error}</p>
                    )}
                  </div>

                  {/* PFAS-Atlas Results */}
                  <div className={`algorithm-result ${molecule.atlas_success ? 'success' : 'error'}`}>
                    <h6>🗺️ PFAS-Atlas</h6>
                    {molecule.atlas_success ? (
                      <div>
                        <p><strong>First Class:</strong> {molecule.atlas_first_class}</p>
                        <p><strong>Second Class:</strong> {molecule.atlas_second_class}</p>
                        <p><small>Time: {(molecule.atlas_time * 1000).toFixed(2)}ms</small></p>
                      </div>
                    ) : (
                      <p className="text-danger">Error: {molecule.atlas_error}</p>
                    )}
                  </div>
                </div>

                {/* Manual Review Panel */}
                <div className="review-panel">
                  <h6>Manual Review</h6>
                  
                  {/* Review Status */}
                  {molecule.review_date && (
                    <Alert variant="info" className="mb-3">
                      Previously reviewed on {new Date(molecule.review_date).toLocaleDateString()}
                    </Alert>
                  )}

                  {/* PFASGroups Correctness */}
                  <Form.Group className="mb-3">
                    <Form.Label>PFASGroups Classification Correct?</Form.Label>
                    <div className="review-buttons">
                      <Button
                        variant={reviewData[molecule.id]?.pfasgroupsCorrect === true ? 'success' : 'outline-success'}
                        onClick={() => handleReviewChange(molecule.id, 'pfasgroupsCorrect', true)}
                      >
                        ✅ Correct
                      </Button>
                      <Button
                        variant={reviewData[molecule.id]?.pfasgroupsCorrect === false ? 'danger' : 'outline-danger'}
                        onClick={() => handleReviewChange(molecule.id, 'pfasgroupsCorrect', false)}
                      >
                        ❌ Incorrect
                      </Button>
                      <Button
                        variant={reviewData[molecule.id]?.pfasgroupsCorrect === null ? 'secondary' : 'outline-secondary'}
                        onClick={() => handleReviewChange(molecule.id, 'pfasgroupsCorrect', null)}
                      >
                        🤷 Unclear
                      </Button>
                    </div>
                  </Form.Group>

                  {/* PFAS-Atlas Correctness */}
                  <Form.Group className="mb-3">
                    <Form.Label>PFAS-Atlas Classification Correct?</Form.Label>
                    <div className="review-buttons">
                      <Button
                        variant={reviewData[molecule.id]?.atlasCorrect === true ? 'success' : 'outline-success'}
                        onClick={() => handleReviewChange(molecule.id, 'atlasCorrect', true)}
                      >
                        ✅ Correct
                      </Button>
                      <Button
                        variant={reviewData[molecule.id]?.atlasCorrect === false ? 'danger' : 'outline-danger'}
                        onClick={() => handleReviewChange(molecule.id, 'atlasCorrect', false)}
                      >
                        ❌ Incorrect
                      </Button>
                      <Button
                        variant={reviewData[molecule.id]?.atlasCorrect === null ? 'secondary' : 'outline-secondary'}
                        onClick={() => handleReviewChange(molecule.id, 'atlasCorrect', null)}
                      >
                        🤷 Unclear
                      </Button>
                    </div>
                  </Form.Group>

                  {/* Notes */}
                  <Form.Group className="mb-3">
                    <Form.Label>Reviewer Notes</Form.Label>
                    <Form.Control
                      as="textarea"
                      rows={3}
                      value={reviewData[molecule.id]?.reviewerNotes || ''}
                      onChange={(e) => handleReviewChange(molecule.id, 'reviewerNotes', e.target.value)}
                      placeholder="Add your notes about this classification..."
                    />
                  </Form.Group>

                  {/* Submit Review */}
                  <div className="d-flex justify-content-between align-items-center">
                    <Button
                      variant="primary"
                      onClick={() => submitReview(molecule.id)}
                      disabled={submitStatus[molecule.id] === 'submitting'}
                    >
                      {submitStatus[molecule.id] === 'submitting' ? (
                        <>
                          <Spinner as="span" animation="border" size="sm" className="me-2" />
                          Submitting...
                        </>
                      ) : (
                        'Submit Review'
                      )}
                    </Button>

                    {submitStatus[molecule.id] === 'success' && (
                      <Alert variant="success" className="mb-0 py-1">Review saved!</Alert>
                    )}
                    {submitStatus[molecule.id] === 'error' && (
                      <Alert variant="danger" className="mb-0 py-1">Error saving review</Alert>
                    )}
                  </div>
                </div>
              </Col>
            </Row>
          </Card.Body>
        </Card>
      ))}

      {/* Pagination */}
      {totalPages > 1 && (
        <div className="pagination-wrapper">
          <ReactPaginate
            pageCount={totalPages}
            onPageChange={handlePageChange}
            containerClassName="pagination"
            activeClassName="active"
            disabledClassName="disabled"
            previousLabel="← Previous"
            nextLabel="Next →"
            forcePage={currentPage - 1}
          />
        </div>
      )}
    </div>
  );
}

export default MoleculeReviewer;