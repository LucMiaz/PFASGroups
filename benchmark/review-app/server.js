const express = require('express');
const cors = require('cors');
const bodyParser = require('body-parser');
const path = require('path');
const fs = require('fs');
const https = require('https');
const Database = require('./database/database');
const ExcelJS = require('exceljs');

const app = express();
const PORT = process.env.PORT || 5000;

// Initialize database
const db = new Database();

// Load PFAS groups mapping
let pfasGroupsMap = {};
try {
    const groupsData = JSON.parse(fs.readFileSync(path.join(__dirname, 'pfas_groups_map.json'), 'utf8'));
    pfasGroupsMap = Object.fromEntries(groupsData.map(g => [g.id, g.name]));
    console.log(`Loaded ${Object.keys(pfasGroupsMap).length} PFAS groups`);
} catch (error) {
    console.error('Warning: Could not load PFAS groups map:', error.message);
}

// Helper function to convert group IDs to names
function enrichGroupData(groupIds) {
    if (!Array.isArray(groupIds)) return [];
    return groupIds.map(id => ({
        id: id,
        name: pfasGroupsMap[id] || `Group ${id}`
    }));
}

// Middleware
app.use(cors());
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({ extended: true }));

// Serve static files from React build
app.use(express.static(path.join(__dirname, 'client/build')));

// Serve analysis reports static files (images, HTML)
app.use('/analysis-reports', express.static(path.join(__dirname, 'analysis_reports')));

// API Routes

// Get molecules with pagination and filtering
app.get('/api/molecules', async (req, res) => {
    try {
        const page = parseInt(req.query.page) || 1;
        const limit = parseInt(req.query.limit) || 20;
        const offset = (page - 1) * limit;
        const dataset = req.query.dataset || 'all';
        const reviewStatus = req.query.reviewStatus || 'all'; // all, reviewed, unreviewed
        const search = req.query.search || '';
        const prioritizeMisclassified = req.query.prioritizeMisclassified !== 'false'; // default true

        // Build WHERE clause
        let whereConditions = [];
        let params = [];

        if (dataset !== 'all') {
            whereConditions.push('m.dataset_type = ?');
            params.push(dataset);
        }

        if (search) {
            whereConditions.push('(m.smiles LIKE ? OR m.group_name LIKE ?)');
            params.push(`%${search}%`, `%${search}%`);
        }

        // Review status filtering
        if (reviewStatus === 'reviewed') {
            whereConditions.push('mr.id IS NOT NULL');
        } else if (reviewStatus === 'unreviewed') {
            whereConditions.push('mr.id IS NULL');
        }

        const whereClause = whereConditions.length > 0 ? `WHERE ${whereConditions.join(' AND ')}` : '';

        // Build ORDER BY clause to prioritize misclassified unreviewed entries and flavor mismatches
        let orderClause = 'ORDER BY m.id ASC';
        if (prioritizeMisclassified) {
            orderClause = `ORDER BY 
                CASE 
                    WHEN mr.id IS NULL AND (pg.detected_groups != pgbc.detected_groups) THEN 1  -- Unreviewed with flavor mismatch (highest priority)
                    WHEN mr.id IS NULL AND (pg.success = 0 OR ar.success = 0) THEN 2            -- Unreviewed misclassified
                    WHEN mr.id IS NULL THEN 3                                                     -- Other unreviewed
                    WHEN pg.detected_groups != pgbc.detected_groups THEN 4                       -- Reviewed with flavor mismatch
                    WHEN mr.pfasgroups_correct = 0 OR mr.atlas_correct = 0 THEN 5                -- Reviewed misclassified
                    ELSE 6                                                                        -- Reviewed correct (lowest priority)
                END,
                m.id ASC`;
        }

        // Get molecules with results including bycomponent flavor
        const molecules = await db.all(`
            SELECT 
                m.*,
                pg.detected_groups as pfasgroups_detected,
                pg.detected_definitions as pfasgroups_detected_definitions,
                pg.success as pfasgroups_success,
                pg.execution_time as pfasgroups_time,
                pg.error_message as pfasgroups_error,
                pgbc.detected_groups as pfasgroups_bycomponent_detected,
                pgbc.detected_definitions as pfasgroups_bycomponent_detected_definitions,
                pgbc.success as pfasgroups_bycomponent_success,
                pgbc.execution_time as pfasgroups_bycomponent_time,
                pgbc.error_message as pfasgroups_bycomponent_error,
                ar.first_class as atlas_first_class,
                ar.second_class as atlas_second_class,
                ar.success as atlas_success,
                ar.execution_time as atlas_time,
                ar.error_message as atlas_error,
                mr.pfasgroups_correct,
                mr.atlas_correct,
                mr.reviewer_notes,
                mr.review_date,
                mr.is_pfas as manual_is_pfas,
                mr.correct_groups as manual_correct_groups,
                mr.correct_classification as manual_correct_classification,
                CASE 
                    WHEN mr.id IS NULL AND (pg.detected_groups != pgbc.detected_groups) THEN 'flavor_mismatch_unreviewed'
                    WHEN mr.id IS NULL AND (pg.success = 0 OR ar.success = 0) THEN 'misclassified_unreviewed'
                    WHEN mr.id IS NULL THEN 'unreviewed'
                    WHEN pg.detected_groups != pgbc.detected_groups THEN 'flavor_mismatch_reviewed'
                    WHEN mr.pfasgroups_correct = 0 OR mr.atlas_correct = 0 THEN 'misclassified_reviewed'
                    ELSE 'reviewed_correct'
                END as priority_category
            FROM molecules m
            LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
            LEFT JOIN pfasgroups_results_bycomponent pgbc ON m.id = pgbc.molecule_id
            LEFT JOIN atlas_results ar ON m.id = ar.molecule_id
            LEFT JOIN manual_reviews mr ON m.id = mr.molecule_id
            ${whereClause}
            ${orderClause}
            LIMIT ? OFFSET ?
        `, [...params, limit, offset]);

        // Get total count for pagination
        const totalCountResult = await db.get(`
            SELECT COUNT(DISTINCT m.id) as total
            FROM molecules m
            LEFT JOIN manual_reviews mr ON m.id = mr.molecule_id
            ${whereClause}
        `, params);

        const totalCount = totalCountResult.total;
        const totalPages = Math.ceil(totalCount / limit);

        // Parse JSON fields and enrich with group names
        const processedMolecules = molecules.map(mol => ({
            ...mol,
            target_groups: mol.target_groups ? JSON.parse(mol.target_groups) : [],
            pfasgroups_detected: enrichGroupData(mol.pfasgroups_detected ? JSON.parse(mol.pfasgroups_detected) : []),
            pfasgroups_detected_definitions: mol.pfasgroups_detected_definitions ? JSON.parse(mol.pfasgroups_detected_definitions) : [],
            pfasgroups_bycomponent_detected: enrichGroupData(mol.pfasgroups_bycomponent_detected ? JSON.parse(mol.pfasgroups_bycomponent_detected) : []),
            pfasgroups_bycomponent_detected_definitions: mol.pfasgroups_bycomponent_detected_definitions ? JSON.parse(mol.pfasgroups_bycomponent_detected_definitions) : [],
            manual_correct_groups: mol.manual_correct_groups ? JSON.parse(mol.manual_correct_groups) : null,
            has_flavor_mismatch: mol.pfasgroups_detected !== mol.pfasgroups_bycomponent_detected
        }));

        res.json({
            molecules: processedMolecules,
            pagination: {
                page,
                limit,
                totalCount,
                totalPages,
                hasNext: page < totalPages,
                hasPrev: page > 1
            }
        });
    } catch (error) {
        console.error('Error fetching molecules:', error);
        res.status(500).json({ error: 'Failed to fetch molecules' });
    }
});

// Get dataset statistics
app.get('/api/stats', async (req, res) => {
    try {
        // Get molecule counts by dataset
        const datasetStats = await db.all(`
            SELECT 
                dataset_type,
                COUNT(*) as total_molecules,
                COUNT(mr.id) as reviewed_molecules
            FROM molecules m
            LEFT JOIN manual_reviews mr ON m.id = mr.molecule_id
            GROUP BY dataset_type
        `);

        // Get overall accuracy stats
        const accuracyStats = await db.all(`
            SELECT
                COUNT(*) as total_reviewed,
                SUM(CASE WHEN pfasgroups_correct = 1 THEN 1 ELSE 0 END) as pfasgroups_correct_count,
                SUM(CASE WHEN atlas_correct = 1 THEN 1 ELSE 0 END) as atlas_correct_count,
                AVG(CASE WHEN pfasgroups_correct IS NOT NULL THEN pfasgroups_correct ELSE NULL END) as pfasgroups_accuracy,
                AVG(CASE WHEN atlas_correct IS NOT NULL THEN atlas_correct ELSE NULL END) as atlas_accuracy
            FROM manual_reviews
        `);

        res.json({
            datasets: datasetStats,
            accuracy: accuracyStats[0] || {}
        });
    } catch (error) {
        console.error('Error fetching stats:', error);
        res.status(500).json({ error: 'Failed to fetch statistics' });
    }
});

// Submit manual review
app.post('/api/review', async (req, res) => {
    try {
        const {
            moleculeId,
            pfasgroupsCorrect,
            atlasCorrect,
            reviewerNotes,
            reviewerName,
            isPfas,
            correctGroups,
            correctClassification
        } = req.body;

        // Check if review already exists
        const existingReview = await db.get(`
            SELECT id FROM manual_reviews WHERE molecule_id = ?
        `, [moleculeId]);

        if (existingReview) {
            // Update existing review
            await db.run(`
                UPDATE manual_reviews 
                SET pfasgroups_correct = ?, atlas_correct = ?, reviewer_notes = ?,
                    reviewer_name = ?, is_pfas = ?, correct_groups = ?,
                    correct_classification = ?, review_date = CURRENT_TIMESTAMP
                WHERE molecule_id = ?
            `, [
                pfasgroupsCorrect,
                atlasCorrect,
                reviewerNotes,
                reviewerName,
                isPfas,
                JSON.stringify(correctGroups || []),
                correctClassification,
                moleculeId
            ]);
        } else {
            // Insert new review
            await db.run(`
                INSERT INTO manual_reviews (
                    molecule_id, pfasgroups_correct, atlas_correct, reviewer_notes,
                    reviewer_name, is_pfas, correct_groups, correct_classification
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            `, [
                moleculeId,
                pfasgroupsCorrect,
                atlasCorrect,
                reviewerNotes,
                reviewerName,
                isPfas,
                JSON.stringify(correctGroups || []),
                correctClassification
            ]);
        }

        res.json({ success: true, message: 'Review saved successfully' });
    } catch (error) {
        console.error('Error saving review:', error);
        res.status(500).json({ error: 'Failed to save review' });
    }
});

// Export reviews to CSV/JSON
app.get('/api/export/reviews', async (req, res) => {
    try {
        const format = req.query.format || 'json';
        
        const reviews = await db.all(`
            SELECT 
                m.id as molecule_id,
                m.smiles,
                m.dataset_type,
                m.group_name,
                m.target_groups,
                pg.detected_groups as pfasgroups_detected,
                pg.detected_definitions as pfasgroups_detected_definitions,
                ar.first_class as atlas_first_class,
                ar.second_class as atlas_second_class,
                mr.pfasgroups_correct,
                mr.atlas_correct,
                mr.reviewer_notes,
                mr.reviewer_name,
                mr.review_date,
                mr.is_pfas,
                mr.correct_groups,
                mr.correct_classification
            FROM manual_reviews mr
            JOIN molecules m ON m.id = mr.molecule_id
            LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
            LEFT JOIN atlas_results ar ON m.id = ar.molecule_id
            ORDER BY mr.review_date DESC
        `);

        const processedReviews = reviews.map(review => ({
            ...review,
            target_groups: review.target_groups ? JSON.parse(review.target_groups) : [],
            pfasgroups_detected: enrichGroupData(review.pfasgroups_detected ? JSON.parse(review.pfasgroups_detected) : []),
            pfasgroups_detected_definitions: review.pfasgroups_detected_definitions ? JSON.parse(review.pfasgroups_detected_definitions) : [],
            correct_groups: review.correct_groups ? JSON.parse(review.correct_groups) : null
        }));

        if (format === 'csv') {
            // Convert to CSV
            const csv = convertToCSV(processedReviews);
            res.setHeader('Content-Type', 'text/csv');
            res.setHeader('Content-Disposition', 'attachment; filename=reviews.csv');
            res.send(csv);
        } else {
            res.json(processedReviews);
        }
    } catch (error) {
        console.error('Error exporting reviews:', error);
        res.status(500).json({ error: 'Failed to export reviews' });
    }
});

// Get accuracy metrics
app.get('/api/accuracy', async (req, res) => {
    try {
        const dataset = req.query.dataset || 'all';
        
        let whereClause = '';
        let params = [];
        
        if (dataset !== 'all') {
            whereClause = 'WHERE m.dataset_type = ?';
            params = [dataset];
        }

        const accuracy = await db.all(`
            SELECT 
                m.dataset_type,
                COUNT(*) as total_reviewed,
                SUM(CASE WHEN mr.pfasgroups_correct = 1 THEN 1 ELSE 0 END) as pfasgroups_correct,
                SUM(CASE WHEN mr.atlas_correct = 1 THEN 1 ELSE 0 END) as atlas_correct,
                AVG(CASE WHEN mr.pfasgroups_correct IS NOT NULL THEN mr.pfasgroups_correct ELSE NULL END) as pfasgroups_accuracy,
                AVG(CASE WHEN mr.atlas_correct IS NOT NULL THEN mr.atlas_correct ELSE NULL END) as atlas_accuracy
            FROM manual_reviews mr
            JOIN molecules m ON m.id = mr.molecule_id
            ${whereClause}
            GROUP BY m.dataset_type
        `, params);

        res.json(accuracy);
    } catch (error) {
        console.error('Error calculating accuracy:', error);
        res.status(500).json({ error: 'Failed to calculate accuracy' });
    }
});

// Get enhanced performance metrics including manual reviews and fallback to original assessments
app.get('/api/performance-metrics', async (req, res) => {
    try {
        const dataset = req.query.dataset || 'all';
        
        let whereClause = '';
        let params = [];
        
        if (dataset !== 'all') {
            whereClause = 'WHERE m.dataset_type = ?';
            params = [dataset];
        }

        // Get comprehensive performance metrics
        const metrics = await db.all(`
            SELECT 
                m.dataset_type,
                COUNT(*) as total_molecules,
                
                -- Manual review metrics (gold standard)
                COUNT(mr.id) as manually_reviewed_count,
                SUM(CASE WHEN mr.pfasgroups_correct = 1 THEN 1 ELSE 0 END) as pfasgroups_manual_correct,
                SUM(CASE WHEN mr.atlas_correct = 1 THEN 1 ELSE 0 END) as atlas_manual_correct,
                AVG(CASE WHEN mr.pfasgroups_correct IS NOT NULL THEN mr.pfasgroups_correct ELSE NULL END) as pfasgroups_manual_accuracy,
                AVG(CASE WHEN mr.atlas_correct IS NOT NULL THEN mr.atlas_correct ELSE NULL END) as atlas_manual_accuracy,
                
                -- Combined metrics (manual review + original assessment fallback)
                SUM(CASE 
                    WHEN mr.pfasgroups_correct IS NOT NULL THEN mr.pfasgroups_correct
                    WHEN pg.success IS NOT NULL THEN pg.success
                    ELSE 0
                END) as pfasgroups_combined_correct,
                SUM(CASE 
                    WHEN mr.atlas_correct IS NOT NULL THEN mr.atlas_correct
                    WHEN ar.success IS NOT NULL THEN ar.success
                    ELSE 0
                END) as atlas_combined_correct,
                
                AVG(CASE 
                    WHEN mr.pfasgroups_correct IS NOT NULL THEN mr.pfasgroups_correct
                    WHEN pg.success IS NOT NULL THEN pg.success
                    ELSE NULL
                END) as pfasgroups_combined_accuracy,
                AVG(CASE 
                    WHEN mr.atlas_correct IS NOT NULL THEN mr.atlas_correct
                    WHEN ar.success IS NOT NULL THEN ar.success
                    ELSE NULL
                END) as atlas_combined_accuracy,
                
                -- Original algorithm success rates
                SUM(CASE WHEN pg.success = 1 THEN 1 ELSE 0 END) as pfasgroups_original_success,
                SUM(CASE WHEN ar.success = 1 THEN 1 ELSE 0 END) as atlas_original_success,
                AVG(CASE WHEN pg.success IS NOT NULL THEN pg.success ELSE NULL END) as pfasgroups_original_accuracy,
                AVG(CASE WHEN ar.success IS NOT NULL THEN ar.success ELSE NULL END) as atlas_original_accuracy,
                
                -- Coverage metrics
                COUNT(pg.id) as pfasgroups_coverage,
                COUNT(ar.id) as atlas_coverage,
                
                -- Misclassification analysis
                SUM(CASE 
                    WHEN mr.id IS NULL AND (pg.success = 0 OR ar.success = 0) THEN 1 
                    ELSE 0 
                END) as unreviewed_misclassified,
                SUM(CASE 
                    WHEN mr.pfasgroups_correct = 0 OR mr.atlas_correct = 0 THEN 1 
                    ELSE 0 
                END) as reviewed_misclassified
                
            FROM molecules m
            LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
            LEFT JOIN atlas_results ar ON m.id = ar.molecule_id
            LEFT JOIN manual_reviews mr ON m.id = mr.molecule_id
            ${whereClause}
            GROUP BY m.dataset_type
        `, params);

        // Calculate additional derived metrics
        const enhancedMetrics = metrics.map(metric => ({
            ...metric,
            manual_review_coverage: metric.total_molecules > 0 ? (metric.manually_reviewed_count / metric.total_molecules) : 0,
            pfasgroups_coverage_percent: metric.total_molecules > 0 ? (metric.pfasgroups_coverage / metric.total_molecules) : 0,
            atlas_coverage_percent: metric.total_molecules > 0 ? (metric.atlas_coverage / metric.total_molecules) : 0,
            priority_review_needed: metric.unreviewed_misclassified,
            confidence_score: metric.manually_reviewed_count > 10 ? 'high' : metric.manually_reviewed_count > 5 ? 'medium' : 'low'
        }));

        res.json(enhancedMetrics);
    } catch (error) {
        console.error('Error calculating performance metrics:', error);
        res.status(500).json({ error: 'Failed to calculate performance metrics' });
    }
});

// Helper function to convert to CSV
function convertToCSV(data) {
    if (!data.length) return '';
    
    const headers = Object.keys(data[0]);
    const csvContent = [
        headers.join(','),
        ...data.map(row => 
            headers.map(header => {
                let cell = row[header];
                if (Array.isArray(cell)) cell = cell.join(';');
                if (typeof cell === 'string') cell = `"${cell.replace(/"/g, '""')}"`;
                return cell;
            }).join(',')
        )
    ].join('\\n');
    
    return csvContent;
}

// Helper function to fetch molecule image from PubChem service
async function fetchMoleculeImage(smiles) {
    return new Promise((resolve, reject) => {
        const encodedSmiles = encodeURIComponent(smiles);
        // Use PubChem's depiction service - more reliable with valid certificates
        const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodedSmiles}/PNG?image_size=300x300`;
        
        https.get(url, (response) => {
            // Follow redirects
            if (response.statusCode === 301 || response.statusCode === 302) {
                https.get(response.headers.location, (redirectResponse) => {
                    if (redirectResponse.statusCode !== 200) {
                        reject(new Error(`Failed to fetch image: ${redirectResponse.statusCode}`));
                        return;
                    }
                    
                    const chunks = [];
                    redirectResponse.on('data', (chunk) => chunks.push(chunk));
                    redirectResponse.on('end', () => {
                        const imageData = Buffer.concat(chunks);
                        resolve(imageData);
                    });
                }).on('error', (error) => {
                    reject(error);
                });
                return;
            }
            
            if (response.statusCode !== 200) {
                reject(new Error(`Failed to fetch image: ${response.statusCode}`));
                return;
            }
            
            const chunks = [];
            response.on('data', (chunk) => chunks.push(chunk));
            response.on('end', () => {
                const imageData = Buffer.concat(chunks);
                resolve(imageData);
            });
        }).on('error', (error) => {
            reject(error);
        });
    });
}

// Analysis Reports API endpoints
app.get('/api/analysis-reports', async (req, res) => {
    try {
        const reportsDir = path.join(__dirname, 'analysis_reports');
        
        // Try to load reports from files
        const reports = {
            timing: null,
            complex: null,
            enhanced: null
        };
        
        const descriptions = {
            timing: '',
            complex: '',
            enhanced: ''
        };
        
        // Try to load timing report
        try {
            const timingPath = path.join(reportsDir, 'timing_analysis.json');
            if (fs.existsSync(timingPath)) {
                reports.timing = JSON.parse(fs.readFileSync(timingPath, 'utf8'));
            }
            const timingDescPath = path.join(reportsDir, 'timing_description.md');
            if (fs.existsSync(timingDescPath)) {
                descriptions.timing = fs.readFileSync(timingDescPath, 'utf8');
            }
        } catch (err) {
            console.error('Error loading timing report:', err);
        }
        
        // Try to load complex report
        try {
            const complexPath = path.join(reportsDir, 'complex_analysis.json');
            if (fs.existsSync(complexPath)) {
                reports.complex = JSON.parse(fs.readFileSync(complexPath, 'utf8'));
            }
            const complexDescPath = path.join(reportsDir, 'complex_description.md');
            if (fs.existsSync(complexDescPath)) {
                descriptions.complex = fs.readFileSync(complexDescPath, 'utf8');
            }
        } catch (err) {
            console.error('Error loading complex report:', err);
        }
        
        // Try to load enhanced report
        try {
            const enhancedPath = path.join(reportsDir, 'enhanced_analysis.json');
            if (fs.existsSync(enhancedPath)) {
                reports.enhanced = JSON.parse(fs.readFileSync(enhancedPath, 'utf8'));
            }
            const enhancedDescPath = path.join(reportsDir, 'enhanced_description.md');
            if (fs.existsSync(enhancedDescPath)) {
                descriptions.enhanced = fs.readFileSync(enhancedDescPath, 'utf8');
            }
        } catch (err) {
            console.error('Error loading enhanced report:', err);
        }
        
        res.json({ reports, descriptions });
    } catch (error) {
        console.error('Error fetching analysis reports:', error);
        res.status(500).json({ error: 'Failed to fetch analysis reports' });
    }
});

app.post('/api/analysis-reports/description', async (req, res) => {
    try {
        const { reportType, description } = req.body;
        const reportsDir = path.join(__dirname, 'analysis_reports');
        
        // Create directory if it doesn't exist
        if (!fs.existsSync(reportsDir)) {
            fs.mkdirSync(reportsDir, { recursive: true });
        }
        
        const descPath = path.join(reportsDir, `${reportType}_description.md`);
        fs.writeFileSync(descPath, description, 'utf8');
        
        res.json({ success: true });
    } catch (error) {
        console.error('Error saving description:', error);
        res.status(500).json({ error: 'Failed to save description' });
    }
});

app.get('/api/analysis-reports/export/:reportType', async (req, res) => {
    try {
        const { reportType } = req.params;
        const reportsDir = path.join(__dirname, 'analysis_reports');
        
        // Load description
        const descPath = path.join(reportsDir, `${reportType}_description.md`);
        let markdown = '';
        if (fs.existsSync(descPath)) {
            markdown = fs.readFileSync(descPath, 'utf8');
        }
        
        // Load report data
        const reportPath = path.join(reportsDir, `${reportType}_analysis.json`);
        if (fs.existsSync(reportPath)) {
            const report = JSON.parse(fs.readFileSync(reportPath, 'utf8'));
            
            // Add figures section
            if (report.figures && report.figures.length > 0) {
                markdown += '\n\n## Figures\n\n';
                report.figures.forEach((fig, idx) => {
                    markdown += `\n### ${fig.title}\n\n`;
                    
                    // Save figure files
                    const figDir = path.join(reportsDir, 'figures');
                    if (!fs.existsSync(figDir)) {
                        fs.mkdirSync(figDir, { recursive: true });
                    }
                    
                    const svgPath = path.join(figDir, `${reportType}_fig${idx+1}.svg`);
                    const pngPath = path.join(figDir, `${reportType}_fig${idx+1}.png`);
                    
                    // Write SVG
                    if (fig.format === 'svg') {
                        fs.writeFileSync(svgPath, fig.data, 'utf8');
                        markdown += `![${fig.title}](figures/${reportType}_fig${idx+1}.svg)\n\n`;
                    } else if (fig.format === 'png') {
                        // Assume base64 encoded
                        const base64Data = fig.data.replace(/^data:image\/png;base64,/, '');
                        fs.writeFileSync(pngPath, base64Data, 'base64');
                        markdown += `![${fig.title}](figures/${reportType}_fig${idx+1}.png)\n\n`;
                    }
                });
            }
        }
        
        res.setHeader('Content-Type', 'text/markdown');
        res.setHeader('Content-Disposition', `attachment; filename="${reportType}_report.md"`);
        res.send(markdown);
    } catch (error) {
        console.error('Error exporting report:', error);
        res.status(500).json({ error: 'Failed to export report' });
    }
});

// Export data to Excel with images
app.get('/api/export/excel', async (req, res) => {
    try {
        const dataset = req.query.dataset || 'all';
        
        // Query all molecules with their results and reviews
        let whereClause = '';
        let params = [];
        
        if (dataset !== 'all') {
            whereClause = 'WHERE m.dataset_type = ?';
            params = [dataset];
        }
        
        const molecules = await db.all(`
            SELECT 
                m.id,
                m.smiles,
                m.molecular_formula,
                m.dataset_type,
                m.group_name,
                m.target_groups,
                pg.detected_groups as pfasgroups_detected,
                pg.detected_definitions as pfasgroups_detected_definitions,
                pg.success as pfasgroups_success,
                pg.execution_time as pfasgroups_time,
                pg.error_message as pfasgroups_error,
                pgbc.detected_groups as pfasgroups_bycomponent_detected,
                pgbc.detected_definitions as pfasgroups_bycomponent_detected_definitions,
                ar.first_class as atlas_first_class,
                ar.second_class as atlas_second_class,
                ar.success as atlas_success,
                ar.execution_time as atlas_time,
                ar.error_message as atlas_error,
                mr.pfasgroups_correct,
                mr.atlas_correct,
                mr.reviewer_notes,
                mr.reviewer_name,
                mr.review_date,
                mr.is_pfas as manual_is_pfas,
                mr.correct_groups as manual_correct_groups,
                mr.correct_classification as manual_correct_classification
            FROM molecules m
            LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
            LEFT JOIN pfasgroups_results_bycomponent pgbc ON m.id = pgbc.molecule_id
            LEFT JOIN atlas_results ar ON m.id = ar.molecule_id
            LEFT JOIN manual_reviews mr ON m.id = mr.molecule_id
            ${whereClause}
            ORDER BY m.id ASC
        `, params);
        
        // Create workbook
        const workbook = new ExcelJS.Workbook();
        workbook.creator = 'PFAS Benchmark Reviewer';
        workbook.created = new Date();
        
        // Add Summary Sheet
        const summarySheet = workbook.addWorksheet('Summary');
        summarySheet.columns = [
            { header: 'Metric', key: 'metric', width: 30 },
            { header: 'Value', key: 'value', width: 20 }
        ];
        
        const totalMolecules = molecules.length;
        const reviewedCount = molecules.filter(m => m.pfasgroups_correct !== null).length;
        const pfasgroupsCorrect = molecules.filter(m => m.pfasgroups_correct === 1).length;
        const atlasCorrect = molecules.filter(m => m.atlas_correct === 1).length;
        const pfasgroupsAccuracy = reviewedCount > 0 ? (pfasgroupsCorrect / reviewedCount * 100).toFixed(2) : 'N/A';
        const atlasAccuracy = reviewedCount > 0 ? (atlasCorrect / reviewedCount * 100).toFixed(2) : 'N/A';
        
        summarySheet.addRows([
            { metric: 'Dataset', value: dataset },
            { metric: 'Total Molecules', value: totalMolecules },
            { metric: 'Reviewed Molecules', value: reviewedCount },
            { metric: 'PFASGroups Correct', value: pfasgroupsCorrect },
            { metric: 'PFASGroups Accuracy', value: pfasgroupsAccuracy + '%' },
            { metric: 'PFAS-Atlas Correct', value: atlasCorrect },
            { metric: 'PFAS-Atlas Accuracy', value: atlasAccuracy + '%' },
            { metric: 'Export Date', value: new Date().toISOString() }
        ]);
        
        summarySheet.getRow(1).font = { bold: true };
        summarySheet.getColumn('metric').font = { bold: true };
        
        // Add Data Sheet with images
        const dataSheet = workbook.addWorksheet('Molecules');
        dataSheet.columns = [
            { header: 'ID', key: 'id', width: 8 },
            { header: 'SMILES', key: 'smiles', width: 50 },
            { header: 'Molecular Formula', key: 'formula', width: 25 },
            { header: 'Structure', key: 'structure', width: 30 },
            { header: 'Dataset', key: 'dataset', width: 20 },
            { header: 'Manual PFAS Classification', key: 'expected_pfas', width: 25 },
            { header: 'PFASGroups Detected', key: 'pfasgroups_detected', width: 40 },
            { header: 'PFASGroups Success', key: 'pfasgroups_success', width: 18 },
            { header: 'PFASGroups Time (ms)', key: 'pfasgroups_time', width: 20 },
            { header: 'Atlas First Class', key: 'atlas_first', width: 20 },
            { header: 'Atlas Second Class', key: 'atlas_second', width: 20 },
            { header: 'Atlas Success', key: 'atlas_success', width: 15 },
            { header: 'Atlas Time (ms)', key: 'atlas_time', width: 18 },
            { header: 'Reviewed', key: 'reviewed', width: 12 },
            { header: 'PFASGroups Correct', key: 'pfasgroups_correct', width: 20 },
            { header: 'Atlas Correct', key: 'atlas_correct', width: 15 },
            { header: 'Reviewer Notes', key: 'notes', width: 40 },
            { header: 'Reviewer', key: 'reviewer', width: 20 },
            { header: 'Review Date', key: 'review_date', width: 20 }
        ];
        
        // Style header row
        dataSheet.getRow(1).font = { bold: true };
        dataSheet.getRow(1).fill = {
            type: 'pattern',
            pattern: 'solid',
            fgColor: { argb: 'FFE0E0E0' }
        };
        
        // Add data rows
        for (let i = 0; i < molecules.length; i++) {
            const mol = molecules[i];
            
            // Parse JSON fields
            const pfasgroupsDetected = mol.pfasgroups_detected ? 
                enrichGroupData(JSON.parse(mol.pfasgroups_detected)).map(g => g.name).join(', ') : '';
            
            const rowData = {
                id: mol.id,
                smiles: mol.smiles,
                formula: mol.molecular_formula || '',
                structure: '', // Placeholder for image
                dataset: mol.dataset_type,
                expected_pfas: mol.manual_is_pfas !== null ? (mol.manual_is_pfas ? 'Yes' : 'No') : 'N/A',
                pfasgroups_detected: pfasgroupsDetected,
                pfasgroups_success: mol.pfasgroups_success ? 'Success' : 'Failed',
                pfasgroups_time: mol.pfasgroups_time ? mol.pfasgroups_time.toFixed(2) : '',
                atlas_first: mol.atlas_first_class || '',
                atlas_second: mol.atlas_second_class || '',
                atlas_success: mol.atlas_success ? 'Success' : 'Failed',
                atlas_time: mol.atlas_time ? mol.atlas_time.toFixed(2) : '',
                reviewed: mol.pfasgroups_correct !== null ? 'Yes' : 'No',
                pfasgroups_correct: mol.pfasgroups_correct === 1 ? 'Correct' : mol.pfasgroups_correct === 0 ? 'Incorrect' : '',
                atlas_correct: mol.atlas_correct === 1 ? 'Correct' : mol.atlas_correct === 0 ? 'Incorrect' : '',
                notes: mol.reviewer_notes || '',
                reviewer: mol.reviewer_name || '',
                review_date: mol.review_date || ''
            };
            
            const row = dataSheet.addRow(rowData);
            
            // Set row height for image
            row.height = 100;
            
            // Try to add molecule image
            try {
                // First try local file
                const imagePath = path.join(__dirname, 'molecule_images', `${mol.id}.png`);
                let imageId;
                
                if (fs.existsSync(imagePath)) {
                    imageId = workbook.addImage({
                        filename: imagePath,
                        extension: 'png',
                    });
                } else if (mol.smiles) {
                    // Fetch image from PubChem service
                    try {
                        const imageData = await fetchMoleculeImage(mol.smiles);
                        imageId = workbook.addImage({
                            buffer: imageData,
                            extension: 'png',
                        });
                    } catch (fetchError) {
                        console.log(`Could not fetch image for molecule ${mol.id}: ${fetchError.message}`);
                    }
                }
                
                if (imageId !== undefined) {
                    dataSheet.addImage(imageId, {
                        tl: { col: 3, row: i + 1 },
                        ext: { width: 180, height: 90 }
                    });
                }
            } catch (imgError) {
                console.log(`Error adding image for molecule ${mol.id}: ${imgError.message}`);
            }
        }
        
        // Add conditional formatting for success/failure
        dataSheet.eachRow((row, rowNumber) => {
            if (rowNumber > 1) { // Skip header
                // PFASGroups Success
                const pgSuccessCell = row.getCell('pfasgroups_success');
                if (pgSuccessCell.value === 'Success') {
                    pgSuccessCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFD4EDDA' }
                    };
                } else if (pgSuccessCell.value === 'Failed') {
                    pgSuccessCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFF8D7DA' }
                    };
                }
                
                // Atlas Success
                const atlasSuccessCell = row.getCell('atlas_success');
                if (atlasSuccessCell.value === 'Success') {
                    atlasSuccessCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFD4EDDA' }
                    };
                } else if (atlasSuccessCell.value === 'Failed') {
                    atlasSuccessCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFF8D7DA' }
                    };
                }
                
                // Review status
                const pgCorrectCell = row.getCell('pfasgroups_correct');
                if (pgCorrectCell.value === 'Correct') {
                    pgCorrectCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFD4EDDA' }
                    };
                    pgCorrectCell.font = { color: { argb: 'FF155724' } };
                } else if (pgCorrectCell.value === 'Incorrect') {
                    pgCorrectCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFF8D7DA' }
                    };
                    pgCorrectCell.font = { color: { argb: 'FF721C24' } };
                }
                
                const atlasCorrectCell = row.getCell('atlas_correct');
                if (atlasCorrectCell.value === 'Correct') {
                    atlasCorrectCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFD4EDDA' }
                    };
                    atlasCorrectCell.font = { color: { argb: 'FF155724' } };
                } else if (atlasCorrectCell.value === 'Incorrect') {
                    atlasCorrectCell.fill = {
                        type: 'pattern',
                        pattern: 'solid',
                        fgColor: { argb: 'FFF8D7DA' }
                    };
                    atlasCorrectCell.font = { color: { argb: 'FF721C24' } };
                }
            }
        });
        
        // Generate Excel file
        const buffer = await workbook.xlsx.writeBuffer();
        
        const filename = `pfas_benchmark_${dataset}_${new Date().toISOString().split('T')[0]}.xlsx`;
        res.setHeader('Content-Type', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet');
        res.setHeader('Content-Disposition', `attachment; filename="${filename}"`);
        res.send(buffer);
        
    } catch (error) {
        console.error('Error exporting to Excel:', error);
        res.status(500).json({ error: 'Failed to export to Excel' });
    }
});

// Catch all handler for React routing
app.get('*', (req, res) => {
    res.sendFile(path.join(__dirname, 'client/build', 'index.html'));
});

// Error handling middleware
app.use((error, req, res, next) => {
    console.error('Server error:', error);
    res.status(500).json({ error: 'Internal server error' });
});

// Start server
app.listen(PORT, () => {
    console.log(`🚀 PFAS Benchmark Reviewer server running on port ${PORT}`);
    console.log(`📊 Database initialized`);
});

// Graceful shutdown
process.on('SIGINT', async () => {
    console.log('\\n🛑 Shutting down server...');
    await db.close();
    process.exit(0);
});

module.exports = app;