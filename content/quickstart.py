from pathlib import Path
import streamlit as st

from src.common.common import page_setup, v_space

page_setup(page="main")


def inject_custom_css():
    """Inject CSS for custom styling with responsive design."""
    st.markdown(
        """
        <style>
        /* Hero section styling */
        .hero-section {
            text-align: center;
            margin-bottom: 2rem;
        }

        .hero-title {
            font-size: 2.5rem;
            font-weight: 700;
            color: #29379b;
            margin-bottom: 0.5rem;
        }

        .hero-subtitle {
            font-size: 1.25rem;
            color: #6c757d;
            margin-bottom: 1.5rem;
        }

        /* Feature card styling */
        .feature-card {
            background: linear-gradient(135deg, #f8f9fa 0%, #f1f3f4 100%);
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            padding: 1.25rem;
            margin: 0.5rem 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }

        .feature-card h4 {
            color: #29379b;
            margin-bottom: 0.5rem;
            font-size: 1.1rem;
        }

        .feature-card p {
            color: #6c757d;
            margin: 0;
            font-size: 0.95rem;
        }

        /* Citation box */
        .citation-box {
            background: linear-gradient(135deg, #f8f9fa 0%, #f1f3f4 100%);
            border-left: 4px solid #29379b;
            border-radius: 0 8px 8px 0;
            padding: 1rem 1.25rem;
            margin: 1rem 0;
        }

        .citation-box p {
            margin: 0;
            color: #495057;
            font-size: 0.95rem;
        }

        /* Main action button styling */
        .st-key-start_workflow_btn button {
            background: linear-gradient(135deg, #29379b 0%, #1e2a7a 100%) !important;
            border: none !important;
            border-radius: 12px !important;
            padding: 1.5rem 3rem !important;
            font-size: 1.25rem !important;
            font-weight: 700 !important;
            color: white !important;
            box-shadow: 0 4px 12px rgba(41, 55, 155, 0.3) !important;
            transition: all 0.3s ease !important;
        }

        .st-key-start_workflow_btn button:hover {
            background: linear-gradient(135deg, #1e2a7a 0%, #162159 100%) !important;
            transform: translateY(-2px) !important;
            box-shadow: 0 6px 16px rgba(41, 55, 155, 0.4) !important;
        }

        .st-key-start_workflow_btn button p {
            color: white !important;
            font-size: 1.25rem !important;
            font-weight: 700 !important;
        }

        /* Download box styling */
        .st-key-download_container {
            background: linear-gradient(135deg, #f8f9fa 0%, #f1f3f4 100%) !important;
            border: 1px solid #e0e0e0 !important;
            border-radius: 8px !important;
            padding: 1.5rem !important;
            margin: 1rem 0 !important;
            text-align: center !important;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05) !important;
        }

        .st-key-download_container > div {
            background: transparent !important;
        }

        /* Responsive design */
        @media (max-width: 768px) {
            .hero-title {
                font-size: 2rem;
            }

            .hero-subtitle {
                font-size: 1.1rem;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_hero_section():
    """Render the hero section with title and logo."""
    st.markdown('<div class="hero-section">', unsafe_allow_html=True)

    title_col, logo_col, spacer2 = st.columns([4, 1.5, 0.5])

    with title_col:
        st.markdown(
            """
            <h1 class="hero-title">MHCquant</h1>
            <p class="hero-subtitle">Scalable and Sensitive HLA Peptide Identification for Immunopeptidomics</p>
            """,
            unsafe_allow_html=True,
        )

    with logo_col:
        st.image("assets/openms_transparent_bg_logo.svg", width=150)

    st.markdown("</div>", unsafe_allow_html=True)


def render_description():
    """Render the main description section."""
    st.markdown(
        """
        MHCquant is an open-source pipeline for **highly sensitive HLA peptide identification** from
        mass spectrometry data. By integrating peptide property predictors **DeepLC** and **MS2PIP**
        via the **MS2Rescore** framework, MHCquant improves peptide identifications by up to **27%**
        compared to conventional approaches, particularly enriching **low-abundant peptides** critical
        for **tumor antigen discovery** and **cancer immunotherapy**.
        """
    )


def render_workflow_overview():
    """Render the workflow steps overview."""
    st.markdown("### Workflow Overview")

    cols = st.columns(5)

    steps = [
        ("1", "Upload", "MS data (mzML) and protein database (FASTA)", "content/topp_workflow_file_upload.py"),
        ("2", "Search", "Database search with Comet for peptide identification", "content/topp_workflow_parameter.py"),
        ("3", "Rescore", "MS2Rescore with DeepLC + MS2PIP, FDR via Percolator", "content/topp_workflow_parameter.py"),
        ("4", "Visualize", "Interactive exploration of identified HLA peptides", "content/topp_workflow_results.py"),
        ("5", "Export", "Download results in TSV and mzTab formats", "content/download_section.py"),
    ]

    for col, (num, title, desc, page) in zip(cols, steps):
        with col:
            button_key = f"workflow_step_{num}_btn"
            if st.button(
                f"{title}",
                key=button_key,
                use_container_width=True,
                help=desc
            ):
                st.switch_page(page)

            # Apply custom styling to make it look like a workflow step card
            st.markdown(
                f"""
                <style>
                .st-key-{button_key} button {{
                    background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%) !important;
                    border: 2px solid #dee2e6 !important;
                    border-radius: 12px !important;
                    padding: 1.5rem 0.5rem !important;
                    min-height: 140px !important;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1) !important;
                    transition: all 0.3s ease !important;
                    display: flex !important;
                    flex-direction: column !important;
                    justify-content: center !important;
                    align-items: center !important;
                }}

                .st-key-{button_key} button::before {{
                    content: "{num}";
                    background: #29379b;
                    color: white;
                    width: 2rem;
                    height: 2rem;
                    border-radius: 50%;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    font-weight: 700;
                    margin-bottom: 0.5rem;
                    font-size: 1rem;
                }}

                .st-key-{button_key} button p {{
                    color: #29379b !important;
                    font-weight: 600 !important;
                    font-size: 1rem !important;
                    margin: 0 !important;
                }}

                .st-key-{button_key} button::after {{
                    content: "{desc}";
                    display: block;
                    font-size: 0.8rem !important;
                    font-weight: 400 !important;
                    color: #6c757d !important;
                    margin-top: 0.5rem !important;
                    text-align: center;
                    line-height: 1.3;
                }}

                .st-key-{button_key} button:hover {{
                    border-color: #29379b !important;
                    transform: translateY(-2px) !important;
                    box-shadow: 0 4px 12px rgba(41, 55, 155, 0.2) !important;
                    background: linear-gradient(135deg, #e9ecef 0%, #dee2e6 100%) !important;
                }}

                .st-key-{button_key} button:hover p {{
                    color: #1e2a7a !important;
                }}
                </style>
                """,
                unsafe_allow_html=True,
            )


def render_features():
    """Render the key features section."""
    st.markdown("### Key Features")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown(
            """
            <div class="feature-card">
                <h4>Up to 27% More Identifications</h4>
                <p>DeepLC and MS2PIP integration via MS2Rescore dramatically boosts identification rates across diverse MS platforms.</p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.markdown(
            """
            <div class="feature-card">
                <h4>Low-Abundant Peptide Detection</h4>
                <p>Enhanced sensitivity for detecting low-expression antigens critical for neoepitope discovery.</p>
            </div>
            """,
            unsafe_allow_html=True,
        )

    with col2:
        st.markdown(
            """
            <div class="feature-card">
                <h4>Tumor Antigen Discovery</h4>
                <p>Identify mutation-derived neoepitopes and tumor-associated antigens as potential immunotherapy targets.</p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.markdown(
            """
            <div class="feature-card">
                <h4>Reproducible Analysis</h4>
                <p>Built on OpenMS with containerized execution for fully reproducible immunopeptidomics workflows.</p>
            </div>
            """,
            unsafe_allow_html=True,
        )


def render_start_button():
    """Render the main start workflow button."""
    st.markdown("---")

    spacer1, center, spacer2 = st.columns([1, 2, 1])

    with center:
        if st.button(
            "Start MHCquant Analysis",
            key="start_workflow_btn",
            use_container_width=True,
            help="Navigate to the file upload page to begin your analysis"
        ):
            st.switch_page("content/topp_workflow_file_upload.py")


def render_download_section():
    """Render the download section for offline app."""
    if not Path("OpenMS-App.zip").exists():
        return

    st.markdown("---")
    st.markdown(
        """
        <h4 style="color: #6c757d; text-align: center; margin-bottom: 1rem;">
            Want to use MHCquant offline?
        </h4>
        """,
        unsafe_allow_html=True,
    )

    spacer1, center, spacer2 = st.columns([2, 2, 2])

    with center:
        with st.container(key="download_container"):
            st.markdown(
                """
                <p style="color: #6c757d; margin-bottom: 1rem;">
                    Download the Windows version for offline analysis.
                </p>
                """,
                unsafe_allow_html=True,
            )

            with open("OpenMS-App.zip", "rb") as file:
                st.download_button(
                    label="Download for Windows",
                    data=file,
                    file_name="MHCquant-App.zip",
                    mime="archive/zip",
                    type="secondary",
                    use_container_width=True,
                )


def render_citation():
    """Render the citation section."""
    st.markdown("---")
    st.markdown("### Citation")
    st.markdown(
        """
        <div class="citation-box">
            <p>
                Scheld JS, Bichmann L, Nelde A, et al. <strong>MHCquant2 refines immunopeptidomics tumor antigen
                discovery.</strong> <em>Genome Biology</em> 2025, 26, 63.
                <a href="https://doi.org/10.1186/s13059-025-03763-8" target="_blank">doi:10.1186/s13059-025-03763-8</a>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="citation-box">
            <p>
                Müller TD, Siraj A, Schmid S, et al. <strong>OpenMS WebApps: Building User-Friendly Solutions
                for MS Analysis.</strong> <em>Journal of Proteome Research</em> 2025.
                <a href="https://doi.org/10.1021/acs.jproteome.4c00872" target="_blank">doi:10.1021/acs.jproteome.4c00872</a>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="citation-box">
            <p>
                Pfeuffer J, Bielow C, Wein S, et al. <strong>OpenMS 3 enables reproducible analysis of large-scale
                mass spectrometry data.</strong> <em>Nature Methods</em> 2024, 21, 365-367.
                <a href="https://doi.org/10.1038/s41592-024-02197-7" target="_blank">doi:10.1038/s41592-024-02197-7</a>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )


def main():
    """Main function to render the quickstart page."""
    inject_custom_css()

    render_hero_section()
    render_description()

    v_space(1)
    render_workflow_overview()

    v_space(1)
    render_features()

    render_start_button()
    render_download_section()
    render_citation()


main()
