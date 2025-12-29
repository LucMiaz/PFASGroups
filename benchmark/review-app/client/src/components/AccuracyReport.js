import React, { useState, useEffect } from 'react';
import { Card, Table, Badge, Button, Row, Col, Alert, ProgressBar, Form } from 'react-bootstrap';
import Select from 'react-select';

function AccuracyReport() {
  const [accuracyData, setAccuracyData] = useState([]);
  const [performanceMetrics, setPerformanceMetrics] = useState([]);
  const [selectedDataset, setSelectedDataset] = useState('all');
  const [loading, setLoading] = useState(true);
  const [exportStatus, setExportStatus] = useState(null);
  const [showEnhancedMetrics, setShowEnhancedMetrics] = useState(true);

  const datasetOptions = [
    { value: 'all', label: 'All Datasets' },
    { value: 'oecd', label: 'OECD' },
    { value: 'enhanced', label: 'Enhanced' },
    { value: 'timing', label: 'Timing' },
    { value: 'complex_branched', label: 'Complex Branched' },
    { value: 'non_fluorinated', label: 'Non-Fluorinated' }
  ];

  useEffect(() => {
    fetchAccuracyData();
  }, [selectedDataset]);

  const fetchAccuracyData = async () => {
    setLoading(true);
    try {
      const params = selectedDataset !== 'all' ? `?dataset=${selectedDataset}` : '';
      
      // Fetch both accuracy and performance metrics
      const [accuracyResponse, performanceResponse] = await Promise.all([
        fetch(`/api/accuracy${params}`),
        fetch(`/api/performance-metrics${params}`)
      ]);
      
      const accuracyData = await accuracyResponse.json();
      const performanceData = await performanceResponse.json();
      
      setAccuracyData(accuracyData);
      setPerformanceMetrics(performanceData);
    } catch (error) {
      console.error('Error fetching accuracy data:', error);
    }
    setLoading(false);
  };

  const exportReviews = async (format = 'json') => {
    setExportStatus('exporting');
    try {
      const response = await fetch(`/api/export/reviews?format=${format}`);
      
      if (response.ok) {
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `pfas_reviews.${format}`;
        a.click();
        window.URL.revokeObjectURL(url);
        setExportStatus('success');
      } else {
        setExportStatus('error');
      }
    } catch (error) {
      console.error('Error exporting reviews:', error);
      setExportStatus('error');
    }

    setTimeout(() => setExportStatus(null), 3000);
  };

  const exportToExcel = async () => {
    setExportStatus('exporting');
    try {
      const params = selectedDataset !== 'all' ? `?dataset=${selectedDataset}` : '';
      const response = await fetch(`/api/export/excel${params}`);
      
      if (response.ok) {
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        const filename = `pfas_benchmark_${selectedDataset}_${new Date().toISOString().split('T')[0]}.xlsx`;
        a.download = filename;
        a.click();
        window.URL.revokeObjectURL(url);
        setExportStatus('success');
      } else {
        setExportStatus('error');
      }
    } catch (error) {
      console.error('Error exporting to Excel:', error);
      setExportStatus('error');
    }

    setTimeout(() => setExportStatus(null), 3000);
  };

  const calculateOverallAccuracy = () => {
    if (!accuracyData.length) return { pfasgroups: 0, atlas: 0, total: 0 };

    const totalReviewed = accuracyData.reduce((sum, d) => sum + d.total_reviewed, 0);
    const pfasgroupsCorrect = accuracyData.reduce((sum, d) => sum + d.pfasgroups_correct, 0);
    const atlasCorrect = accuracyData.reduce((sum, d) => sum + d.atlas_correct, 0);

    return {
      pfasgroups: totalReviewed > 0 ? (pfasgroupsCorrect / totalReviewed) * 100 : 0,
      atlas: totalReviewed > 0 ? (atlasCorrect / totalReviewed) * 100 : 0,
      total: totalReviewed
    };
  };

  const overall = calculateOverallAccuracy();

  if (loading) {
    return <div>Loading accuracy report...</div>;
  }

  return (
    <div>
      <h1>📈 Accuracy Report</h1>

      {/* Controls */}
      <Card className="mb-4">
        <Card.Body>
          <Row>
            <Col md={3}>
              <h6>Dataset Filter</h6>
              <Select
                value={datasetOptions.find(opt => opt.value === selectedDataset)}
                onChange={(option) => setSelectedDataset(option.value)}
                options={datasetOptions}
              />
            </Col>
            <Col md={2}>
              <h6>View Options</h6>
              <Form.Check
                type="switch"
                id="enhanced-metrics-switch"
                label="Show Enhanced Performance Metrics"
                checked={showEnhancedMetrics}
                onChange={(e) => setShowEnhancedMetrics(e.target.checked)}
              />
            </Col>
            <Col md={3}>
              <h6>Export Reviews</h6>
              <div className="d-flex gap-2">
                <Button 
                  variant="outline-primary" 
                  size="sm"
                  onClick={() => exportReviews('json')}
                  disabled={exportStatus === 'exporting'}
                >
                  Export JSON
                </Button>
                <Button 
                  variant="outline-secondary" 
                  size="sm"
                  onClick={() => exportReviews('csv')}
                  disabled={exportStatus === 'exporting'}
                >
                  Export CSV
                </Button>
              </div>
              {exportStatus === 'success' && (
                <Alert variant="success" className="mt-2 mb-0 py-1">
                  Export completed successfully!
                </Alert>
              )}
              {exportStatus === 'error' && (
                <Alert variant="danger" className="mt-2 mb-0 py-1">
                  Export failed. Please try again.
                </Alert>
              )}
            </Col>
            <Col md={4}>
              <h6>Export to Excel</h6>
              <Button 
                variant="success" 
                size="sm"
                onClick={exportToExcel}
                disabled={exportStatus === 'exporting'}
                className="w-100"
              >
                📊 Export Complete Dataset to Excel
              </Button>
              <small className="text-muted d-block mt-1">
                Includes structure images, formulas, and all classification data
              </small>
            </Col>
          </Row>
        </Card.Body>
      </Card>

      {/* Overall Accuracy Summary */}
      {overall.total > 0 && (
        <Row className="mb-4">
          <Col md={4}>
            <Card className="stats-card">
              <Card.Body className="text-center">
                <h3>{overall.total}</h3>
                <p>Total Reviewed</p>
              </Card.Body>
            </Card>
          </Col>
          <Col md={4}>
            <Card className="stats-card">
              <Card.Body className="text-center">
                <h3>{overall.pfasgroups.toFixed(1)}%</h3>
                <p>PFASGroups Accuracy</p>
              </Card.Body>
            </Card>
          </Col>
          <Col md={4}>
            <Card className="stats-card">
              <Card.Body className="text-center">
                <h3>{overall.atlas.toFixed(1)}%</h3>
                <p>PFAS-Atlas Accuracy</p>
              </Card.Body>
            </Card>
          </Col>
        </Row>
      )}

      {/* Enhanced Performance Metrics */}
      {showEnhancedMetrics && performanceMetrics.length > 0 && (
        <Card className="mb-4">
          <Card.Body>
            <div className="d-flex justify-content-between align-items-center mb-3">
              <h5>Enhanced Performance Metrics</h5>
              <Button 
                variant="outline-secondary" 
                size="sm"
                onClick={() => setShowEnhancedMetrics(!showEnhancedMetrics)}
              >
                Hide Enhanced Metrics
              </Button>
            </div>
            
            <Table responsive striped>
              <thead>
                <tr>
                  <th>Dataset</th>
                  <th>Total Molecules</th>
                  <th>Manual Review Coverage</th>
                  <th>Manual Accuracy</th>
                  <th>Combined Accuracy</th>
                  <th>Priority Reviews Needed</th>
                  <th>Confidence</th>
                </tr>
              </thead>
              <tbody>
                {performanceMetrics.map((metric, index) => (
                  <tr key={index}>
                    <td>
                      <Badge bg="info">{metric.dataset_type.toUpperCase()}</Badge>
                    </td>
                    <td>{metric.total_molecules}</td>
                    <td>
                      <div className="d-flex align-items-center">
                        <div className="me-2" style={{width: '100px'}}>
                          <ProgressBar 
                            now={metric.manual_review_coverage * 100} 
                            variant={metric.manual_review_coverage > 0.8 ? 'success' : metric.manual_review_coverage > 0.5 ? 'warning' : 'danger'}
                            size="sm"
                          />
                        </div>
                        <small>{(metric.manual_review_coverage * 100).toFixed(1)}% ({metric.manually_reviewed_count}/{metric.total_molecules})</small>
                      </div>
                    </td>
                    <td>
                      <div>
                        <Badge bg="primary" className="me-1">
                          PG: {metric.pfasgroups_manual_accuracy ? (metric.pfasgroups_manual_accuracy * 100).toFixed(1) : 'N/A'}%
                        </Badge>
                        <br/>
                        <Badge bg="info">
                          PA: {metric.atlas_manual_accuracy ? (metric.atlas_manual_accuracy * 100).toFixed(1) : 'N/A'}%
                        </Badge>
                      </div>
                    </td>
                    <td>
                      <div>
                        <Badge bg="success" className="me-1">
                          PG: {metric.pfasgroups_combined_accuracy ? (metric.pfasgroups_combined_accuracy * 100).toFixed(1) : 'N/A'}%
                        </Badge>
                        <br/>
                        <Badge bg="warning">
                          PA: {metric.atlas_combined_accuracy ? (metric.atlas_combined_accuracy * 100).toFixed(1) : 'N/A'}%
                        </Badge>
                      </div>
                    </td>
                    <td>
                      <Badge bg={metric.priority_review_needed > 0 ? 'danger' : 'success'}>
                        {metric.priority_review_needed} molecules
                      </Badge>
                    </td>
                    <td>
                      <Badge bg={metric.confidence_score === 'high' ? 'success' : metric.confidence_score === 'medium' ? 'warning' : 'danger'}>
                        {metric.confidence_score.toUpperCase()}
                      </Badge>
                    </td>
                  </tr>
                ))}
              </tbody>
            </Table>
            
            <Alert variant="info" className="mt-3">
              <small>
                <strong>Legend:</strong><br/>
                • <strong>Manual Accuracy:</strong> Based only on manually reviewed molecules (gold standard)<br/>
                • <strong>Combined Accuracy:</strong> Manual review when available, falls back to original algorithm assessment<br/>
                • <strong>PG:</strong> PFASGroups, <strong>PA:</strong> PFAS-Atlas<br/>
                • <strong>Priority Reviews:</strong> Unreviewed molecules that were misclassified by algorithms
              </small>
            </Alert>
          </Card.Body>
        </Card>
      )}

      {/* Detailed Accuracy Table */}
      <Card>
        <Card.Body>
          <h5>Accuracy by Dataset</h5>
          
          {accuracyData.length === 0 ? (
            <Alert variant="info">
              No reviewed molecules found for the selected dataset. 
              Start reviewing molecules to see accuracy metrics.
            </Alert>
          ) : (
            <Table responsive striped>
              <thead>
                <tr>
                  <th>Dataset</th>
                  <th>Total Reviewed</th>
                  <th>PFASGroups</th>
                  <th>PFAS-Atlas</th>
                  <th>Comparison</th>
                </tr>
              </thead>
              <tbody>
                {accuracyData.map((dataset, index) => {
                  const pfasgroupsAcc = dataset.pfasgroups_accuracy ? (dataset.pfasgroups_accuracy * 100) : 0;
                  const atlasAcc = dataset.atlas_accuracy ? (dataset.atlas_accuracy * 100) : 0;
                  
                  return (
                    <tr key={index}>
                      <td>
                        <Badge bg="info">{dataset.dataset_type.toUpperCase()}</Badge>
                      </td>
                      <td>{dataset.total_reviewed}</td>
                      <td>
                        <div className="d-flex align-items-center">
                          <div className="accuracy-badge">
                            <Badge bg={pfasgroupsAcc >= 90 ? 'success' : pfasgroupsAcc >= 70 ? 'warning' : 'danger'}>
                              {pfasgroupsAcc.toFixed(1)}%
                            </Badge>
                          </div>
                          <small className="ms-2 text-muted">
                            ({dataset.pfasgroups_correct}/{dataset.total_reviewed})
                          </small>
                        </div>
                      </td>
                      <td>
                        <div className="d-flex align-items-center">
                          <div className="accuracy-badge">
                            <Badge bg={atlasAcc >= 90 ? 'success' : atlasAcc >= 70 ? 'warning' : 'danger'}>
                              {atlasAcc.toFixed(1)}%
                            </Badge>
                          </div>
                          <small className="ms-2 text-muted">
                            ({dataset.atlas_correct}/{dataset.total_reviewed})
                          </small>
                        </div>
                      </td>
                      <td>
                        {pfasgroupsAcc > atlasAcc ? (
                          <Badge bg="primary">PFASGroups +{(pfasgroupsAcc - atlasAcc).toFixed(1)}%</Badge>
                        ) : atlasAcc > pfasgroupsAcc ? (
                          <Badge bg="info">PFAS-Atlas +{(atlasAcc - pfasgroupsAcc).toFixed(1)}%</Badge>
                        ) : (
                          <Badge bg="secondary">Equal</Badge>
                        )}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </Table>
          )}
        </Card.Body>
      </Card>
    </div>
  );
}

export default AccuracyReport;