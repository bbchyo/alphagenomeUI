"""
🧬 AlphaGenome DNA Analiz Aracı / DNA Analysis Tool

Modern Streamlit arayüzü ile DNA sekans analizi, varyant efekt tahmini ve genomik bölge analizi.

Powered by Google DeepMind AlphaGenome API
https://deepmind.google.com/science/alphagenome

Geliştirici / Developer: Arş. Gör. Hakkı Halil BABACAN
Bilimin gelişimini hızlandırmak için açık kaynak olarak sunulmuştur.

⚠️ YASAL UYARI / LEGAL NOTICE:
- Bu araç yalnızca akademik ve araştırma amaçlı kullanım içindir.
- Ticari kullanım yasaktır.
- Tüm çıktılar AlphaGenome Output Terms of Use şartlarına tabidir:
  https://deepmind.google.com/science/alphagenome/output-terms

⚠️ NOT FOR COMMERCIAL USE - Academic and research purposes only.
All outputs are subject to AlphaGenome Output Terms of Use.
"""

import streamlit as st
import json
import base64
from datetime import datetime

# AlphaGenome SDK imports (conditionally imported)
try:
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    ALPHAGENOME_AVAILABLE = True
except ImportError:
    ALPHAGENOME_AVAILABLE = False

# =============================================================================
# PAGE CONFIGURATION
# =============================================================================
st.set_page_config(
    page_title="AlphaGenome - DNA Analiz Aracı",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# =============================================================================
# TRANSLATIONS
# =============================================================================
TRANSLATIONS = {
    "tr": {
        "title": "🧬 AlphaGenome DNA Analiz Aracı",
        "subtitle": "Google DeepMind AlphaGenome modeli ile genomik analizler",
        "api_key": "API Anahtarı",
        "api_key_placeholder": "AlphaGenome API anahtarınızı girin",
        "api_key_help": "API anahtarınızı [AlphaGenome Portal](https://deepmind.google.com/science/alphagenome) adresinden alabilirsiniz.",
        "language": "Dil",
        "sequence_analysis": "Sekans Analizi",
        "variant_analysis": "Varyant Tahmini",
        "interval_analysis": "Bölge Analizi",
        "sequence": "DNA Sekansı",
        "sequence_placeholder": "DNA sekansınızı girin (A, T, G, C, N)...",
        "sequence_length": "Sekans Uzunluğu",
        "min_length_warning": "⚠️ Minimum 16,384 baz çifti gerekli",
        "organism": "Organizma",
        "human": "İnsan",
        "mouse": "Fare",
        "tissue": "Doku/Hücre Tipi",
        "output_type": "Çıktı Türü",
        "chromosome": "Kromozom",
        "position": "Pozisyon",
        "reference": "Referans Baz",
        "alternate": "Alternatif Baz",
        "start": "Başlangıç",
        "end": "Bitiş",
        "analyze": "Analiz Et",
        "load_example": "Örnek Yükle",
        "clear": "Temizle",
        "results": "Sonuçlar",
        "download_json": "JSON İndir",
        "analyzing": "Analiz yapılıyor...",
        "error": "Hata",
        "success": "Başarılı",
        "no_api_key": "Lütfen önce API anahtarınızı girin",
        "developer": "Geliştirici",
        "contribute": "Katkıda Bulunun",
        "contribute_text": "Bilimin gelişimini hızlandırmak için bu projeye katkıda bulunabilirsiniz. Kollektif bilime inanıyoruz!",
        "non_commercial": "⚠️ Ticari amaçlar için kullanılamaz. Yalnızca akademik ve araştırma amaçlıdır.",
        "sdk_not_installed": "⚠️ AlphaGenome SDK yüklü değil. Lütfen `pip install git+https://github.com/google-deepmind/alphagenome.git` komutu ile yükleyin.",
    },
    "en": {
        "title": "🧬 AlphaGenome DNA Analysis Tool",
        "subtitle": "Genomic analysis with Google DeepMind AlphaGenome model",
        "api_key": "API Key",
        "api_key_placeholder": "Enter your AlphaGenome API key",
        "api_key_help": "Get your API key from [AlphaGenome Portal](https://deepmind.google.com/science/alphagenome).",
        "language": "Language",
        "sequence_analysis": "Sequence Analysis",
        "variant_analysis": "Variant Prediction",
        "interval_analysis": "Region Analysis",
        "sequence": "DNA Sequence",
        "sequence_placeholder": "Enter your DNA sequence (A, T, G, C, N)...",
        "sequence_length": "Sequence Length",
        "min_length_warning": "⚠️ Minimum 16,384 base pairs required",
        "organism": "Organism",
        "human": "Human",
        "mouse": "Mouse",
        "tissue": "Tissue/Cell Type",
        "output_type": "Output Type",
        "chromosome": "Chromosome",
        "position": "Position",
        "reference": "Reference Base",
        "alternate": "Alternate Base",
        "start": "Start",
        "end": "End",
        "analyze": "Analyze",
        "load_example": "Load Example",
        "clear": "Clear",
        "results": "Results",
        "download_json": "Download JSON",
        "analyzing": "Analyzing...",
        "error": "Error",
        "success": "Success",
        "no_api_key": "Please enter your API key first",
        "developer": "Developer",
        "contribute": "Contribute",
        "contribute_text": "Contribute to this project to accelerate scientific progress. We believe in collective science!",
        "non_commercial": "⚠️ Not for commercial use. For academic and research purposes only.",
        "sdk_not_installed": "⚠️ AlphaGenome SDK not installed. Please run `pip install git+https://github.com/google-deepmind/alphagenome.git`",
    }
}

TISSUES = {
    "UBERON:0002048": {"tr": "Akciğer", "en": "Lung"},
    "UBERON:0000955": {"tr": "Beyin", "en": "Brain"},
    "UBERON:0002107": {"tr": "Karaciğer", "en": "Liver"},
    "UBERON:0002106": {"tr": "Dalak", "en": "Spleen"},
    "UBERON:0002113": {"tr": "Böbrek", "en": "Kidney"},
    "UBERON:0000948": {"tr": "Kalp", "en": "Heart"},
    "UBERON:0001157": {"tr": "Kolon", "en": "Colon"},
}

OUTPUT_TYPES = {
    "RNA_SEQ": {"tr": "RNA-seq (Gen Ekspresyonu)", "en": "RNA-seq (Gene Expression)"},
    "DNASE": {"tr": "DNase-seq (DNA Erişilebilirliği)", "en": "DNase-seq (DNA Accessibility)"},
    "ATAC": {"tr": "ATAC-seq (Kromatin Erişilebilirliği)", "en": "ATAC-seq (Chromatin Accessibility)"},
    "CAGE": {"tr": "CAGE-seq (Transkripsiyon Başlangıcı)", "en": "CAGE-seq (Transcription Start)"},
    "CHIP_HISTONE": {"tr": "ChIP-seq Histon Modifikasyonları", "en": "ChIP-seq Histone Modifications"},
}

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
BASES = ["A", "T", "G", "C"]

# Example sequence for testing
EXAMPLE_SEQUENCE = "A" * 8192 + "TGCA" * 2048

# =============================================================================
# CUSTOM CSS
# =============================================================================
def load_css():
    st.markdown("""
    <style>
    /* Dark theme */
    .stApp {
        background: linear-gradient(135deg, #0f0f23 0%, #1a1a3e 50%, #0f172a 100%);
    }
    
    /* Header */
    .main-header {
        background: linear-gradient(135deg, rgba(99, 102, 241, 0.1), rgba(16, 185, 129, 0.1));
        border: 1px solid rgba(99, 102, 241, 0.3);
        border-radius: 16px;
        padding: 2rem;
        margin-bottom: 2rem;
        text-align: center;
    }
    
    .main-header h1 {
        background: linear-gradient(135deg, #818cf8, #34d399);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 2.5rem;
        margin-bottom: 0.5rem;
    }
    
    /* Cards */
    .glass-card {
        background: rgba(30, 41, 59, 0.8);
        backdrop-filter: blur(16px);
        border: 1px solid rgba(99, 102, 241, 0.2);
        border-radius: 16px;
        padding: 1.5rem;
        margin-bottom: 1rem;
    }
    
    /* Footer */
    .footer {
        background: rgba(15, 23, 42, 0.9);
        border-top: 1px solid rgba(99, 102, 241, 0.2);
        padding: 2rem;
        margin-top: 3rem;
        text-align: center;
        border-radius: 16px;
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #6366f1, #10b981);
        color: white;
        border: none;
        border-radius: 12px;
        padding: 0.75rem 2rem;
        font-weight: 600;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 25px rgba(99, 102, 241, 0.4);
    }
    
    /* Inputs */
    .stTextInput > div > div > input,
    .stTextArea > div > div > textarea,
    .stSelectbox > div > div > select {
        background: rgba(30, 41, 59, 0.8);
        border: 1px solid rgba(99, 102, 241, 0.3);
        border-radius: 8px;
        color: white;
    }
    
    /* Success/Error messages */
    .success-box {
        background: rgba(16, 185, 129, 0.2);
        border: 1px solid rgba(16, 185, 129, 0.5);
        border-radius: 12px;
        padding: 1rem;
        margin: 1rem 0;
    }
    
    .error-box {
        background: rgba(239, 68, 68, 0.2);
        border: 1px solid rgba(239, 68, 68, 0.5);
        border-radius: 12px;
        padding: 1rem;
        margin: 1rem 0;
    }
    
    /* DNA animation */
    @keyframes pulse {
        0%, 100% { opacity: 0.5; }
        50% { opacity: 1; }
    }
    
    .dna-icon {
        animation: pulse 2s ease-in-out infinite;
    }
    </style>
    """, unsafe_allow_html=True)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def get_translation(key: str, lang: str = "tr") -> str:
    """Get translation for a key."""
    return TRANSLATIONS.get(lang, TRANSLATIONS["tr"]).get(key, key)

def t(key: str) -> str:
    """Shorthand for getting translation with current language."""
    return get_translation(key, st.session_state.get("language", "tr"))

def validate_sequence(sequence: str) -> bool:
    """Validate DNA sequence contains only valid characters."""
    valid_chars = set("ATGCN")
    return all(c.upper() in valid_chars for c in sequence.replace(" ", "").replace("\n", ""))

def get_output_type_enum(output_type: str):
    """Get the AlphaGenome OutputType enum value."""
    if ALPHAGENOME_AVAILABLE:
        type_map = {
            "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
            "DNASE": dna_client.OutputType.DNASE,
            "ATAC": dna_client.OutputType.ATAC,
            "CAGE": dna_client.OutputType.CAGE,
            "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
        }
        return type_map.get(output_type, dna_client.OutputType.RNA_SEQ)
    return None

def create_download_link(data: dict, filename: str) -> str:
    """Create a download link for JSON data."""
    json_str = json.dumps(data, indent=2, default=str)
    b64 = base64.b64encode(json_str.encode()).decode()
    return f'<a href="data:application/json;base64,{b64}" download="{filename}" class="download-btn">📥 {t("download_json")}</a>'

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================
def analyze_sequence(api_key: str, sequence: str, organism: str, tissue: str, output_type: str):
    """Perform DNA sequence analysis."""
    try:
        model = dna_client.create(api_key)
        
        outputs = model.predict_sequence(
            sequence=sequence,
            organism=organism,
            ontology_terms=[tissue],
            requested_outputs=[get_output_type_enum(output_type)],
        )
        
        result = {
            "type": "sequence",
            "sequence_length": len(sequence),
            "organism": organism,
            "tissue": tissue,
            "output_type": output_type,
            "data_shape": str(getattr(outputs, output_type.lower(), outputs).values.shape) if hasattr(outputs, output_type.lower()) else "N/A",
            "timestamp": datetime.now().isoformat(),
        }
        
        return {"success": True, "data": result}
    except Exception as e:
        return {"success": False, "error": str(e)}

def analyze_variant(api_key: str, chromosome: str, position: int, ref: str, alt: str, tissue: str, output_type: str):
    """Perform variant effect prediction."""
    try:
        model = dna_client.create(api_key)
        
        # Create interval around variant (need context)
        interval = genome.Interval(
            chromosome=chromosome,
            start=max(0, position - 500000),
            end=position + 500000
        )
        
        variant = genome.Variant(
            chromosome=chromosome,
            position=position,
            reference_bases=ref,
            alternate_bases=alt,
        )
        
        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=[tissue],
            requested_outputs=[get_output_type_enum(output_type)],
        )
        
        result = {
            "type": "variant",
            "chromosome": chromosome,
            "position": position,
            "reference": ref,
            "alternate": alt,
            "tissue": tissue,
            "output_type": output_type,
            "ref_shape": str(outputs.reference.rna_seq.values.shape) if hasattr(outputs.reference, 'rna_seq') else "N/A",
            "alt_shape": str(outputs.alternate.rna_seq.values.shape) if hasattr(outputs.alternate, 'rna_seq') else "N/A",
            "timestamp": datetime.now().isoformat(),
        }
        
        return {"success": True, "data": result}
    except Exception as e:
        return {"success": False, "error": str(e)}

def analyze_interval(api_key: str, chromosome: str, start: int, end: int, tissue: str, output_type: str):
    """Perform genomic interval analysis."""
    try:
        model = dna_client.create(api_key)
        
        interval = genome.Interval(
            chromosome=chromosome,
            start=start,
            end=end
        )
        
        outputs = model.predict_interval(
            interval=interval,
            ontology_terms=[tissue],
            requested_outputs=[get_output_type_enum(output_type)],
        )
        
        result = {
            "type": "interval",
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "tissue": tissue,
            "output_type": output_type,
            "data_shape": str(outputs.rna_seq.values.shape) if hasattr(outputs, 'rna_seq') else "N/A",
            "timestamp": datetime.now().isoformat(),
        }
        
        return {"success": True, "data": result}
    except Exception as e:
        return {"success": False, "error": str(e)}

# =============================================================================
# UI COMPONENTS
# =============================================================================
def render_sidebar():
    """Render the sidebar with settings."""
    with st.sidebar:
        st.markdown("---")
        
        # Language selector
        lang = st.selectbox(
            t("language"),
            options=["tr", "en"],
            format_func=lambda x: "🇹🇷 Türkçe" if x == "tr" else "🇬🇧 English",
            key="language"
        )
        
        st.markdown("---")
        
        # API Key
        st.markdown(f"### 🔑 {t('api_key')}")
        api_key = st.text_input(
            t("api_key"),
            type="password",
            placeholder=t("api_key_placeholder"),
            key="api_key",
            label_visibility="collapsed"
        )
        st.markdown(t("api_key_help"))
        
        st.markdown("---")
        
        # Developer info
        st.markdown(f"### 👨‍💻 {t('developer')}")
        st.markdown("""
        **Arş. Gör. Hakkı Halil BABACAN**
        """)
        
        st.markdown("---")
        
        # Contribute
        st.markdown(f"### 🤝 {t('contribute')}")
        st.markdown(t("contribute_text"))
        st.markdown("🔗 [GitHub](https://github.com/bbchyo/alphagenomeUI)")
        
        # Non-commercial warning
        st.warning(t("non_commercial"))
        
        return api_key

def render_header():
    """Render the main header."""
    st.markdown(f"""
    <div class="main-header">
        <h1>{t("title")}</h1>
        <p style="color: #94a3b8; font-size: 1.1rem;">{t("subtitle")}</p>
    </div>
    """, unsafe_allow_html=True)

def render_sequence_tab(api_key: str):
    """Render the sequence analysis tab."""
    col1, col2 = st.columns([2, 1])
    
    # Initialize sequence in session state if needed
    if "seq_text" not in st.session_state:
        st.session_state["seq_text"] = ""
    
    with col1:
        sequence = st.text_area(
            t("sequence"),
            value=st.session_state.get("seq_text", ""),
            placeholder=t("sequence_placeholder"),
            height=200,
            key="sequence_input_area"
        )
        
        # Update session state
        st.session_state["seq_text"] = sequence
        
        # Show sequence length
        if sequence:
            clean_seq = sequence.replace(" ", "").replace("\n", "")
            length = len(clean_seq)
            if length < 16384:
                st.warning(f"{t('sequence_length')}: {length:,} bp - {t('min_length_warning')}")
            else:
                st.success(f"{t('sequence_length')}: {length:,} bp ✓")
    
    with col2:
        organism = st.selectbox(
            t("organism"),
            options=["human", "mouse"],
            format_func=lambda x: t("human") if x == "human" else t("mouse")
        )
        
        tissue_options = list(TISSUES.keys())
        tissue = st.selectbox(
            t("tissue"),
            options=tissue_options,
            format_func=lambda x: TISSUES[x][st.session_state.get("language", "tr")]
        )
        
        output_options = list(OUTPUT_TYPES.keys())
        output_type = st.selectbox(
            t("output_type"),
            options=output_options,
            format_func=lambda x: OUTPUT_TYPES[x][st.session_state.get("language", "tr")]
        )
        
        col_btn1, col_btn2 = st.columns(2)
        with col_btn1:
            if st.button(t("load_example"), key="seq_example"):
                st.session_state["seq_text"] = EXAMPLE_SEQUENCE
                st.rerun()
        with col_btn2:
            if st.button(t("clear"), key="seq_clear"):
                st.session_state["seq_text"] = ""
                st.rerun()
    
    if st.button(t("analyze"), key="analyze_sequence", type="primary", use_container_width=True):
        if not api_key:
            st.error(t("no_api_key"))
        elif not sequence or len(sequence.replace(" ", "").replace("\n", "")) < 16384:
            st.error(t("min_length_warning"))
        elif not validate_sequence(sequence):
            st.error("Invalid DNA sequence")
        else:
            if ALPHAGENOME_AVAILABLE:
                with st.spinner(t("analyzing")):
                    result = analyze_sequence(api_key, sequence.replace(" ", "").replace("\n", ""), organism, tissue, output_type)
                    if result["success"]:
                        st.success(f"✅ {t('success')}")
                        st.json(result["data"])
                        st.markdown(create_download_link(result["data"], "sequence_analysis.json"), unsafe_allow_html=True)
                    else:
                        st.error(f"❌ {t('error')}: {result['error']}")
            else:
                st.error(t("sdk_not_installed"))

def render_variant_tab(api_key: str):
    """Render the variant analysis tab."""
    col1, col2, col3 = st.columns(3)
    
    with col1:
        chromosome = st.selectbox(t("chromosome"), options=CHROMOSOMES)
        position = st.number_input(t("position"), min_value=1, value=36201698)
    
    with col2:
        ref = st.selectbox(t("reference"), options=BASES)
        alt = st.selectbox(t("alternate"), options=BASES, index=1)
    
    with col3:
        tissue_options = list(TISSUES.keys())
        tissue = st.selectbox(
            t("tissue"),
            options=tissue_options,
            format_func=lambda x: TISSUES[x][st.session_state.get("language", "tr")],
            key="variant_tissue"
        )
        
        output_options = list(OUTPUT_TYPES.keys())
        output_type = st.selectbox(
            t("output_type"),
            options=output_options,
            format_func=lambda x: OUTPUT_TYPES[x][st.session_state.get("language", "tr")],
            key="variant_output"
        )
    
    if st.button(t("load_example"), key="variant_example"):
        st.info("Example: chr22:36201698 A → C (UBERON:0001157 - Colon)")
    
    if st.button(t("analyze"), key="analyze_variant", type="primary", use_container_width=True):
        if not api_key:
            st.error(t("no_api_key"))
        elif ref == alt:
            st.error("Reference and alternate bases must be different")
        else:
            if ALPHAGENOME_AVAILABLE:
                with st.spinner(t("analyzing")):
                    result = analyze_variant(api_key, chromosome, position, ref, alt, tissue, output_type)
                    if result["success"]:
                        st.success(f"✅ {t('success')}")
                        st.json(result["data"])
                        st.markdown(create_download_link(result["data"], "variant_analysis.json"), unsafe_allow_html=True)
                    else:
                        st.error(f"❌ {t('error')}: {result['error']}")
            else:
                st.error(t("sdk_not_installed"))

def render_interval_tab(api_key: str):
    """Render the interval analysis tab."""
    col1, col2 = st.columns(2)
    
    with col1:
        chromosome = st.selectbox(t("chromosome"), options=CHROMOSOMES, key="interval_chr")
        start = st.number_input(t("start"), min_value=0, value=35677410)
        end = st.number_input(t("end"), min_value=1, value=36725986)
    
    with col2:
        tissue_options = list(TISSUES.keys())
        tissue = st.selectbox(
            t("tissue"),
            options=tissue_options,
            format_func=lambda x: TISSUES[x][st.session_state.get("language", "tr")],
            key="interval_tissue"
        )
        
        output_options = list(OUTPUT_TYPES.keys())
        output_type = st.selectbox(
            t("output_type"),
            options=output_options,
            format_func=lambda x: OUTPUT_TYPES[x][st.session_state.get("language", "tr")],
            key="interval_output"
        )
        
        if st.button(t("load_example"), key="interval_example"):
            st.info("Example: chr22:35677410-36725986 (UBERON:0001157 - Colon)")
    
    if st.button(t("analyze"), key="analyze_interval", type="primary", use_container_width=True):
        if not api_key:
            st.error(t("no_api_key"))
        elif start >= end:
            st.error("Start position must be less than end position")
        else:
            if ALPHAGENOME_AVAILABLE:
                with st.spinner(t("analyzing")):
                    result = analyze_interval(api_key, chromosome, start, end, tissue, output_type)
                    if result["success"]:
                        st.success(f"✅ {t('success')}")
                        st.json(result["data"])
                        st.markdown(create_download_link(result["data"], "interval_analysis.json"), unsafe_allow_html=True)
                    else:
                        st.error(f"❌ {t('error')}: {result['error']}") 
            else:
                st.error(t("sdk_not_installed"))

def render_footer():
    """Render the footer."""
    st.markdown("""
    <div class="footer">
        <p style="color: #64748b; font-size: 0.9rem;">Made with ❤️ by Arş. Gör. Hakkı Halil BABACAN</p>
        <p style="color: #475569; font-size: 0.8rem;">Bilimin gelişimini hızlandırmak için açık kaynak olarak sunulmuştur.</p>
        <p style="color: #4285F4; font-size: 0.8rem; margin-top: 0.5rem;">Powered by <a href="https://deepmind.google.com/science/alphagenome" style="color: #4285F4;">Google DeepMind AlphaGenome</a></p>
        <p style="color: #64748b; font-size: 0.7rem; margin-top: 0.5rem;">⚠️ Çıktılar <a href="https://deepmind.google.com/science/alphagenome/output-terms" style="color: #94a3b8;">AlphaGenome Output Terms of Use</a> şartlarına tabidir.</p>
    </div>
    """, unsafe_allow_html=True)

# =============================================================================
# MAIN APP
# =============================================================================
def main():
    """Main application entry point."""
    # Load custom CSS
    load_css()
    
    # Initialize session state
    if "language" not in st.session_state:
        st.session_state["language"] = "tr"
    
    # Render sidebar and get API key
    api_key = render_sidebar()
    
    # Render header
    render_header()
    
    # Check if SDK is available
    if not ALPHAGENOME_AVAILABLE:
        st.error(t("sdk_not_installed"))
        st.code("pip install git+https://github.com/google-deepmind/alphagenome.git", language="bash")
        st.info("After installing, restart the Streamlit app.")
    
    # Create tabs
    tab1, tab2, tab3 = st.tabs([
        f"🧬 {t('sequence_analysis')}",
        f"🔬 {t('variant_analysis')}",
        f"📊 {t('interval_analysis')}"
    ])
    
    with tab1:
        render_sequence_tab(api_key)
    
    with tab2:
        render_variant_tab(api_key)
    
    with tab3:
        render_interval_tab(api_key)
    
    # Render footer
    render_footer()

if __name__ == "__main__":
    main()
