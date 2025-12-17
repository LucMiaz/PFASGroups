#!/usr/bin/env python3

import json
import sys
from generate_unified_report import (
    find_latest_benchmark_files,
    load_benchmark_data,
    create_interactive_time_performance_plot,
    create_interactive_sankey_diagram
)

try:
    from enhanced_analysis import create_enhanced_sankey_comparison, analyze_system_comparison
    enhanced_analysis_available = True
except ImportError:
    enhanced_analysis_available = False

def test_interactive_plots():
    print("🔍 Testing interactive plot generation...")
    
    files = find_latest_benchmark_files()
    print(f"📂 Found files: {files}")
    
    # Test timing plot
    if files['timing']:
        print("\n⚡ Testing timing plot...")
        timing_data = load_benchmark_data(files['timing'])
        if timing_data:
            timing_html = create_interactive_time_performance_plot(timing_data)
            if timing_html:
                print("✅ Timing plot generated successfully")
                print(f"   HTML length: {len(timing_html)} characters")
                if "plotly" in timing_html.lower():
                    print("✅ Contains Plotly references")
                else:
                    print("❌ No Plotly references found")
            else:
                print("❌ Timing plot generation failed")
        else:
            print("❌ No timing data loaded")
    
    # Test enhanced Sankey diagrams
    if files['enhanced'] and enhanced_analysis_available:
        print("\n🔗 Testing enhanced Sankey diagrams...")
        enhanced_data = load_benchmark_data(files['enhanced'])
        if enhanced_data:
            try:
                # Analyze the enhanced data
                single_analysis, multi_analysis = analyze_system_comparison(enhanced_data)
                # Get the three Sankey diagrams
                sankey_figures = create_enhanced_sankey_comparison(single_analysis, multi_analysis, enhanced_data)
                print(f"   📊 Enhanced Sankey figures returned: {len(sankey_figures) if sankey_figures else 0}")
                
                if sankey_figures and len(sankey_figures) >= 3:
                    for i, fig in enumerate(sankey_figures[:3]):
                        if fig:
                            try:
                                # Generate interactive HTML with unique div IDs
                                include_js = 'cdn' if i == 0 else False
                                sankey_html = fig.to_html(include_plotlyjs=include_js, div_id=f"sankey-{i}")
                                print(f"   ✅ Sankey {i+1} generated: {len(sankey_html)} characters")
                            except Exception as e:
                                print(f"   ❌ Sankey {i+1} HTML conversion failed: {e}")
                        else:
                            print(f"   ❌ Sankey {i+1} figure is None")
                else:
                    print("   ❌ Enhanced Sankey figures not generated properly")
            except Exception as e:
                print(f"   ❌ Enhanced Sankey creation failed: {e}")
                import traceback
                traceback.print_exc()
    
    # Test fallback Sankey diagram
    if files['enhanced']:
        print("\n🔗 Testing fallback Sankey diagram...")
        enhanced_data = load_benchmark_data(files['enhanced'])
        if enhanced_data:
            sankey_html = create_interactive_sankey_diagram(enhanced_data)
            if sankey_html:
                print("✅ Fallback Sankey diagram generated successfully")
                print(f"   HTML length: {len(sankey_html)} characters")
                if "plotly" in sankey_html.lower():
                    print("✅ Contains Plotly references")
                else:
                    print("❌ No Plotly references found")
            else:
                print("❌ Fallback Sankey diagram generation failed")
        else:
            print("❌ No enhanced data loaded")
    
    print("\n🎯 Interactive plot test complete!")

if __name__ == "__main__":
    test_interactive_plots()