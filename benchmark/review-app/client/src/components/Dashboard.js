import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Badge, ProgressBar, Table } from 'react-bootstrap';

function Dashboard({ stats }) {
  const [detailedStats, setDetailedStats] = useState(null);

  useEffect(() => {
    if (stats) {
      calculateDetailedStats();
    }
  }, [stats]);

  const calculateDetailedStats = () => {
    const datasets = stats.datasets || [];
    const accuracy = stats.accuracy || {};

    const totalMolecules = datasets.reduce((sum, d) => sum + d.total_molecules, 0);
    const totalReviewed = datasets.reduce((sum, d) => sum + d.reviewed_molecules, 0);
    const reviewProgress = totalMolecules > 0 ? (totalReviewed / totalMolecules) * 100 : 0;

    setDetailedStats({
      totalMolecules,
      totalReviewed,
      reviewProgress,
      pfasgroupsAccuracy: accuracy.pfasgroups_accuracy ? (accuracy.pfasgroups_accuracy * 100).toFixed(1) : 'N/A',
      atlasAccuracy: accuracy.atlas_accuracy ? (accuracy.atlas_accuracy * 100).toFixed(1) : 'N/A'
    });
  };

  if (!stats || !detailedStats) {
    return <div>Loading dashboard...</div>;
  }

  return (
    <div>
      <h1>📊 PFAS Benchmark Dashboard</h1>
      
      {/* Overview Cards */}
      <Row className="mb-4">
        <Col md={3}>
          <Card className="stats-card">
            <Card.Body>
              <h3>{detailedStats.totalMolecules}</h3>
              <p>Total Molecules</p>
            </Card.Body>
          </Card>
        </Col>
        <Col md={3}>
          <Card className="stats-card">
            <Card.Body>
              <h3>{detailedStats.totalReviewed}</h3>
              <p>Reviewed</p>
            </Card.Body>
          </Card>
        </Col>
        <Col md={3}>
          <Card className="stats-card">
            <Card.Body>
              <h3>{detailedStats.pfasgroupsAccuracy}%</h3>
              <p>PFASGroups Accuracy</p>
            </Card.Body>
          </Card>
        </Col>
        <Col md={3}>
          <Card className="stats-card">
            <Card.Body>
              <h3>{detailedStats.atlasAccuracy}%</h3>
              <p>PFAS-Atlas Accuracy</p>
            </Card.Body>
          </Card>
        </Col>
      </Row>

      {/* Review Progress */}
      <Card className="mb-4">
        <Card.Body>
          <h5>Review Progress</h5>
          <ProgressBar 
            now={detailedStats.reviewProgress} 
            label={`${detailedStats.reviewProgress.toFixed(1)}%`}
            variant={detailedStats.reviewProgress < 30 ? 'danger' : detailedStats.reviewProgress < 70 ? 'warning' : 'success'}
          />
          <small className="text-muted">
            {detailedStats.totalReviewed} of {detailedStats.totalMolecules} molecules reviewed
          </small>
        </Card.Body>
      </Card>

      {/* Dataset Breakdown */}
      <Card>
        <Card.Body>
          <h5>Dataset Breakdown</h5>
          <Table responsive striped>
            <thead>
              <tr>
                <th>Dataset</th>
                <th>Total Molecules</th>
                <th>Reviewed</th>
                <th>Progress</th>
                <th>Status</th>
              </tr>
            </thead>
            <tbody>
              {stats.datasets.map((dataset, index) => {
                const progress = dataset.total_molecules > 0 ? 
                  (dataset.reviewed_molecules / dataset.total_molecules) * 100 : 0;
                
                let variant = 'secondary';
                if (progress > 80) variant = 'success';
                else if (progress > 50) variant = 'warning';
                else if (progress > 0) variant = 'info';

                return (
                  <tr key={index}>
                    <td>
                      <strong>{dataset.dataset_type.toUpperCase()}</strong>
                    </td>
                    <td>{dataset.total_molecules}</td>
                    <td>{dataset.reviewed_molecules}</td>
                    <td>
                      <ProgressBar 
                        now={progress} 
                        variant={variant}
                        style={{height: '20px'}}
                      />
                      <small>{progress.toFixed(1)}%</small>
                    </td>
                    <td>
                      <Badge bg={variant}>
                        {progress === 0 ? 'Not Started' : 
                         progress === 100 ? 'Complete' : 'In Progress'}
                      </Badge>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </Table>
        </Card.Body>
      </Card>
    </div>
  );
}

export default Dashboard;