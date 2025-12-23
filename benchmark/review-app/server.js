const express = require('express');
const cors = require('cors');
const bodyParser = require('body-parser');
const path = require('path');
const Database = require('./database/database');

const app = express();
const PORT = process.env.PORT || 5000;

// Initialize database
const db = new Database();

// Middleware
app.use(cors());
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({ extended: true }));

// Serve static files from React build
app.use(express.static(path.join(__dirname, 'client/build')));

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
                pg.success as pfasgroups_success,
                pg.execution_time as pfasgroups_time,
                pg.error_message as pfasgroups_error,
                pgbc.detected_groups as pfasgroups_bycomponent_detected,
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

        // Parse JSON fields
        const processedMolecules = molecules.map(mol => ({
            ...mol,
            target_groups: mol.target_groups ? JSON.parse(mol.target_groups) : [],
            pfasgroups_detected: mol.pfasgroups_detected ? JSON.parse(mol.pfasgroups_detected) : [],
            pfasgroups_bycomponent_detected: mol.pfasgroups_bycomponent_detected ? JSON.parse(mol.pfasgroups_bycomponent_detected) : [],
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
            pfasgroups_detected: review.pfasgroups_detected ? JSON.parse(review.pfasgroups_detected) : [],
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