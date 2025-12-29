import React, { useState, useEffect } from 'react';
import { BrowserRouter as Router, Routes, Route, Link } from 'react-router-dom';
import { Container, Navbar, Nav } from 'react-bootstrap';
import MoleculeReviewer from './components/MoleculeReviewer';
import Dashboard from './components/Dashboard';
import AccuracyReport from './components/AccuracyReport';
import AnalysisReports from './components/AnalysisReports';
import 'bootstrap/dist/css/bootstrap.min.css';
import './App.css';

function App() {
  const [stats, setStats] = useState(null);

  useEffect(() => {
    fetchStats();
  }, []);

  const fetchStats = async () => {
    try {
      const response = await fetch('/api/stats');
      const data = await response.json();
      setStats(data);
    } catch (error) {
      console.error('Error fetching stats:', error);
    }
  };

  return (
    <Router>
      <div className="App">
        <Navbar bg="dark" variant="dark" expand="lg">
          <Container>
            <Navbar.Brand as={Link} to="/">
              🧪 PFAS Benchmark Reviewer
            </Navbar.Brand>
            <Navbar.Toggle aria-controls="basic-navbar-nav" />
            <Navbar.Collapse id="basic-navbar-nav">
              <Nav className="me-auto">
                <Nav.Link as={Link} to="/">Dashboard</Nav.Link>
                <Nav.Link as={Link} to="/review">Review Molecules</Nav.Link>
                <Nav.Link as={Link} to="/accuracy">Accuracy Report</Nav.Link>
                <Nav.Link as={Link} to="/analysis">Analysis Reports</Nav.Link>
              </Nav>
              {stats && (
                <Navbar.Text>
                  Total Reviewed: {stats.accuracy.total_reviewed || 0}
                </Navbar.Text>
              )}
            </Navbar.Collapse>
          </Container>
        </Navbar>

        <Container className="mt-4">
          <Routes>
            <Route path="/" element={<Dashboard stats={stats} />} />
            <Route path="/review" element={<MoleculeReviewer onReviewUpdate={fetchStats} />} />
            <Route path="/accuracy" element={<AccuracyReport />} />
            <Route path="/analysis" element={<AnalysisReports />} />
          </Routes>
        </Container>
      </div>
    </Router>
  );
}

export default App;