import React, { useState, useEffect } from 'react';
import { Container, Card, Row, Col, Button, Form, Alert, Tabs, Tab } from 'react-bootstrap';
import axios from 'axios';

function AnalysisReports() {
  const [reports, setReports] = useState({
    timing: null,
    complex: null,
    enhanced: null
  });
  const [descriptions, setDescriptions] = useState({
    timing: '',
    complex: '',
    enhanced: ''
  });
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    loadReports();
  }, []);

  const loadReports = async () => {
    setLoading(true);
    try {
      const response = await axios.get('/api/analysis-reports');
      setReports(response.data.reports || {});
      setDescriptions(response.data.descriptions || {
        timing: '',
        complex: '',
        enhanced: ''
      });
      setError(null);
    } catch (err) {
      console.error('Error loading reports:', err);
      setError('Failed to load analysis reports');
    }
    setLoading(false);
  };

  const saveDescription = async (reportType) => {
    try {
      await axios.post('/api/analysis-reports/description', {
        reportType,
        description: descriptions[reportType]
      });
      alert('Description saved successfully!');
    } catch (err) {
      console.error('Error saving description:', err);
      alert('Failed to save description');
    }
  };

  const exportToMarkdown = async (reportType) => {
    try {
      const response = await axios.get(`/api/analysis-reports/export/${reportType}`, {
        responseType: 'blob'
      });
      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', `${reportType}_report.md`);
      document.body.appendChild(link);
      link.click();
      link.remove();
    } catch (err) {
      console.error('Error exporting report:', err);
      alert('Failed to export report');
    }
  };

  const descriptionTemplate = {
    timing: `# Timing Performance Analysis

## Overview
[Describe the overall timing performance characteristics]

## Key Findings
- Finding 1: [Description]
- Finding 2: [Description]
- Finding 3: [Description]

## Performance Characteristics
### PFASGroups
[Describe PFASGroups performance]

### PFAS-Atlas
[Describe PFAS-Atlas performance]

## Scalability Analysis
[Discuss how the algorithms scale with molecule size/complexity]

## Conclusions
[Summary and conclusions]`,

    complex: `# Complex Branched PFAS Analysis

## Overview
[Describe the analysis of complex branched PFAS molecules]

## Test Molecules
[Describe the types of complex molecules tested]

## Detection Results
### PFASGroups
[Describe detection accuracy and patterns]

### PFAS-Atlas  
[Describe detection accuracy and patterns]

## Challenges Identified
- Challenge 1: [Description]
- Challenge 2: [Description]

## Conclusions
[Summary and conclusions]`,

    enhanced: `# Enhanced Functional Groups Analysis

## Overview
[Describe the comprehensive functional group testing]

## Coverage Analysis
[Discuss the breadth of functional groups tested]

## Accuracy Results
### PFASGroups
- Overall accuracy: [X%]
- Strengths: [Description]
- Limitations: [Description]

### PFAS-Atlas
- Overall accuracy: [X%]
- Strengths: [Description]
- Limitations: [Description]

## Comparison
[Compare the two methods]

## Conclusions
[Summary and conclusions]`
  };

  const insertTemplate = (reportType) => {
    setDescriptions({
      ...descriptions,
      [reportType]: descriptionTemplate[reportType]
    });
  };

  const renderReport = (reportType, reportData) => {
    if (!reportData) {
      return (
        <Alert variant="info">
          No {reportType} analysis report available. Run the corresponding Python script to generate it:
          <br />
          <code>python {reportType === 'timing' ? 'analyze_timing.py' : reportType === 'complex' ? 'analyze_complex.py' : 'enhanced_analysis.py'}</code>
        </Alert>
      );
    }

    return (
      <div>
        <Card className="mb-3">
          <Card.Header>
            <strong>Description</strong>
            <Button 
              variant="outline-primary" 
              size="sm" 
              className="float-end ms-2"
              onClick={() => saveDescription(reportType)}
            >
              💾 Save
            </Button>
            <Button 
              variant="outline-secondary" 
              size="sm" 
              className="float-end"
              onClick={() => insertTemplate(reportType)}
            >
              📝 Insert Template
            </Button>
          </Card.Header>
          <Card.Body>
            <Form.Control
              as="textarea"
              rows={10}
              value={descriptions[reportType]}
              onChange={(e) => setDescriptions({
                ...descriptions,
                [reportType]: e.target.value
              })}
              placeholder="Add your analysis description here (Markdown supported)..."
            />
          </Card.Body>
        </Card>

        <Card className="mb-3">
          <Card.Header>
            <strong>Summary Statistics</strong>
          </Card.Header>
          <Card.Body>
            <pre>{JSON.stringify(reportData.summary || reportData, null, 2)}</pre>
          </Card.Body>
        </Card>

        {reportData.figures && reportData.figures.length > 0 && (
          <Card className="mb-3">
            <Card.Header>
              <strong>Figures</strong>
            </Card.Header>
            <Card.Body>
              <Row>
                {reportData.figures.map((figure, idx) => (
                  <Col md={6} key={idx} className="mb-3">
                    <Card>
                      <Card.Header>{figure.title}</Card.Header>
                      <Card.Body>
                        {figure.format === 'svg' ? (
                          <div dangerouslySetInnerHTML={{ __html: figure.data }} />
                        ) : (
                          <img src={figure.data} alt={figure.title} className="img-fluid" />
                        )}
                      </Card.Body>
                    </Card>
                  </Col>
                ))}
              </Row>
            </Card.Body>
          </Card>
        )}

        <div className="text-end">
          <Button 
            variant="success" 
            onClick={() => exportToMarkdown(reportType)}
          >
            📄 Export to Markdown (with figures)
          </Button>
        </div>
      </div>
    );
  };

  if (loading) {
    return (
      <Container className="mt-4">
        <div className="text-center">
          <div className="spinner-border" role="status">
            <span className="visually-hidden">Loading...</span>
          </div>
          <p>Loading analysis reports...</p>
        </div>
      </Container>
    );
  }

  return (
    <Container className="mt-4">
      <h2>📊 Analysis Reports</h2>
      <p className="text-muted">
        View, describe, and export analysis results from timing, complex branched, and enhanced benchmarks.
      </p>

      {error && <Alert variant="danger">{error}</Alert>}

      <Tabs defaultActiveKey="timing" className="mb-3">
        <Tab eventKey="timing" title="⏱️ Timing Performance">
          {renderReport('timing', reports.timing)}
        </Tab>
        <Tab eventKey="complex" title="🧬 Complex Branched">
          {renderReport('complex', reports.complex)}
        </Tab>
        <Tab eventKey="enhanced" title="🔬 Enhanced Analysis">
          {renderReport('enhanced', reports.enhanced)}
        </Tab>
      </Tabs>
    </Container>
  );
}

export default AnalysisReports;
