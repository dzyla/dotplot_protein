import streamlit as st
import numpy as np
import plotly.graph_objects as go
from Bio import SeqIO
from Bio.Align import substitution_matrices
import io
import pandas as pd
from numba import njit

# Set page configuration
st.set_page_config(page_title="Sequence Dot Plot Analysis", layout="wide")

# Initialize session state variables
if 'sequences' not in st.session_state:
    st.session_state.sequences = {'seq1': '', 'seq2': ''}
if 'original_scores' not in st.session_state:
    st.session_state.original_scores = None
if 'filtered_scores' not in st.session_state:
    st.session_state.filtered_scores = None
if 'current_settings' not in st.session_state:
    st.session_state.current_settings = {}
if 'threshold' not in st.session_state:
    st.session_state.threshold = 0.0
if 'window_size' not in st.session_state:
    st.session_state.window_size = 5

# Define DNA scoring matrices as NumPy arrays for faster computation
DNA_MATRIX = {
    'Identity': np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]),
    'Mismatch': np.array([
        [0, -1, -1, -1],
        [-1, 0, -1, -1],
        [-1, -1, 0, -1],
        [-1, -1, -1, 0]
    ])
}

# Mapping for DNA bases to indices
DNA_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}  # 'N' mapped to 'A' as default

# Manually define available protein substitution matrices
PROTEIN_MATRICES = [
    'BLOSUM62',
    'BLOSUM80',
    'PAM250',
    'PAM120',
    'BLOSUM45',
    'BLOSUM50',
    'BLOSUM30',
    'PAM30',
    'PAM70',
    'PAM100'
]

@st.cache_resource
def load_protein_matrix(matrix_name):
    """Load protein substitution matrix using Biopython."""
    try:
        matrix = substitution_matrices.load(matrix_name)
        return matrix
    except ValueError:
        st.error(f"Matrix '{matrix_name}' not found in Biopython substitution matrices.")
        return None

def load_sequence_from_file(file_content, file_name):
    """Load sequence from uploaded file."""
    content = file_content.decode('utf-8').strip()
    if file_name.endswith(('.fasta', '.fa', '.fsa')):
        try:
            records = list(SeqIO.parse(io.StringIO(content), 'fasta'))
            return str(records[0].seq).upper() if records else ''
        except Exception as e:
            st.error(f"Error parsing FASTA file: {e}")
            return ''
    return ''.join(content.split()).upper()

def validate_sequence(seq, seq_type):
    """Validate the sequence based on its type."""
    if seq_type == 'DNA':
        valid_chars = set('ACGTN')
    else:
        # Standard amino acids plus some ambiguity codes
        valid_chars = set('ACDEFGHIKLMNPQRSTVWYBXZJ')
    invalid = set(seq) - valid_chars
    if invalid:
        st.warning(f"Sequence contains invalid characters: {', '.join(sorted(invalid))}")
    return ''.join([c for c in seq if c in valid_chars])

@njit
def calculate_window_scores_dna(seq1_encoded, seq2_encoded, window_size, score_matrix):
    """Optimized window score calculation for DNA using Numba."""
    len1, len2 = len(seq1_encoded), len(seq2_encoded)
    scores = np.zeros((len1 - window_size + 1, len2 - window_size + 1))
    for i in range(len1 - window_size + 1):
        for j in range(len2 - window_size + 1):
            score = 0
            for k in range(window_size):
                score += score_matrix[seq1_encoded[i + k], seq2_encoded[j + k]]
            scores[i, j] = score
    return scores

def calculate_dot_plot(seq1, seq2, window_size, matrix_type, matrix):
    """Calculate dot plot scores using sliding window."""
    if not seq1 or not seq2:
        return None

    len1, len2 = len(seq1), len(seq2)
    if window_size > len1 or window_size > len2:
        st.error("Window size is larger than one of the sequences.")
        return None

    if matrix_type == 'DNA':
        # Encode sequences
        seq1_encoded = np.array([DNA_MAP.get(base, 0) for base in seq1], dtype=np.int32)
        seq2_encoded = np.array([DNA_MAP.get(base, 0) for base in seq2], dtype=np.int32)
        # Calculate scores using Numba-optimized function
        scores = calculate_window_scores_dna(seq1_encoded, seq2_encoded, window_size, DNA_MATRIX[matrix])
    else:
        # For protein sequences, use vectorized operations
        scores = np.zeros((len1 - window_size + 1, len2 - window_size + 1))
        for i in range(len1 - window_size + 1):
            for j in range(len2 - window_size + 1):
                window_score = 0
                for k in range(window_size):
                    a, b = seq1[i + k], seq2[j + k]
                    pair = (a, b)
                    rev_pair = (b, a)
                    window_score += matrix.get(pair, matrix.get(rev_pair, 0))
                scores[i, j] = window_score

    return scores

def apply_threshold(scores, threshold):
    """Apply threshold to scores by setting values below the threshold to 0."""
    return np.where(scores < threshold, 0, scores)

def create_plot(filtered_scores, window_size, threshold, color_scale):
    """Create Plotly heatmap figure."""
    fig = go.Figure(data=go.Heatmap(
        z=filtered_scores,
        colorscale=color_scale,
        colorbar=dict(title='Score'),
        hovertemplate='Score: %{z}<br>Seq1 pos: %{y}<br>Seq2 pos: %{x}<extra></extra>'
    ))
    fig.update_layout(
        title=f'Dot Plot (Window Size: {window_size}, Threshold: {threshold})',
        xaxis_title='Sequence 2 Position',
        yaxis_title='Sequence 1 Position',
        template='plotly_dark',
        autosize=True,
        width=800,
        height=800
    )
    fig.layout.yaxis.scaleanchor = "x"  # Keep aspect ratio square
    return fig

def main():
    st.title('🔍 Sequence Dot Plot Analysis')

    # Sidebar for matrix selection and settings
    st.sidebar.header("🔧 Settings")

    matrix_type = st.sidebar.selectbox(
        'Select Sequence Type:',
        ['DNA', 'Protein'],
        help="Choose the type of sequences you're analyzing."
    )

    if matrix_type == 'DNA':
        matrix_options = list(DNA_MATRIX.keys())
    else:
        matrix_options = PROTEIN_MATRICES

    selected_matrix = st.sidebar.selectbox(
        'Select Scoring Matrix:',
        matrix_options,
        help="Choose a substitution matrix for scoring."
    )

    if matrix_type == 'Protein':
        matrix = load_protein_matrix(selected_matrix)
        if matrix is None:
            st.stop()
    else:
        matrix = selected_matrix  # For DNA, matrix is a string key

    # Main form for sequence input and parameters
    with st.form("sequence_form"):
        st.subheader("📥 Input Sequences")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Sequence 1**")
            seq1_file = st.file_uploader(
                "Upload Sequence 1 (FASTA or TXT)",
                type=['fasta', 'fa', 'fsa', 'txt'],
                key='seq1_file',
                help="Upload a FASTA or plain text file for Sequence 1."
            )
            seq1_text = st.text_area(
                "Or Paste Sequence 1:",
                height=150,
                help="Enter the first sequence here.",
                placeholder='''AACAVDAGSVDQTVQLGQVRTASLAQEGATSSAVGFNIQLNDCDTNVASKAAVAFLGTAIDAGHTNVLALQSSAAGSATNVGVQILDRTGAALTLDGATFSSETTLNNGTNTIPFQARYFATGAATPGAANADATFKVQYQGGGGGGANVVEGKFHVTGGNVTTAA'''
            )

        with col2:
            st.markdown("**Sequence 2**")
            seq2_file = st.file_uploader(
                "Upload Sequence 2 (FASTA or TXT)",
                type=['fasta', 'fa', 'fsa', 'txt'],
                key='seq2_file',
                help="Upload a FASTA or plain text file for Sequence 2.",
            )
            seq2_text = st.text_area(
                "Or Paste Sequence 2:",
                height=150,
                help="Enter the second sequence here.",
                placeholder='''AACAVDAGSVDQTVQLGQVRTASLAQEGATSSAVGFNIQLNDCDTNVASKAAVAFLGTAIDAGHTNVLALQSSAAGSATNVGVQILDRTGAALTLDGATFSSETTLNNGTNTIPFQARYFATGAATPGAANADATFKVQYQGGGGGGANVVEGKFHVTGGNVTTAA'''
            )

        c1, c2, c3 = st.columns(3)
        
        window_size = c1.number_input(
            "🪟 Window Size:",
            min_value=1,
            max_value=50,
            value=st.session_state.window_size,
            step=1,
            help="Size of the sliding window for comparison."
        )

        color_scale = c2.selectbox(
            "🎨 Color Scale:",
            ['Viridis', 'Cividis', 'Plasma', 'Inferno', 'Magma'],
            index=0,
            help="Choose a color scale for the heatmap visualization."
        )

        # Add threshold input within the form
        threshold = c3.number_input(
            "🔽 Threshold Value:",
            min_value=0.0,
            max_value=1000.0,
            value=st.session_state.threshold,
            step=1.0,
            help="Set a threshold value. Scores below this value will be set to 0."
        )

        submitted = st.form_submit_button("Generate Dot Plot")

    if submitted:
        # Load and process sequences
        if seq1_file:
            seq1 = load_sequence_from_file(seq1_file.getvalue(), seq1_file.name)
        else:
            seq1 = seq1_text.strip().upper()

        if seq2_file:
            seq2 = load_sequence_from_file(seq2_file.getvalue(), seq2_file.name)
        else:
            seq2 = seq2_text.strip().upper()

        seq1 = validate_sequence(seq1, matrix_type)
        seq2 = validate_sequence(seq2, matrix_type)

        if not seq1 or not seq2:
            st.error("⚠️ Please provide valid sequences for both Sequence 1 and Sequence 2.")
        else:
            new_settings = {
                'seq1': seq1,
                'seq2': seq2,
                'window_size': window_size,
                'matrix_type': matrix_type,
                'selected_matrix': selected_matrix,
                'color_scale': color_scale
            }

            # Determine if main settings have changed
            if new_settings != st.session_state.current_settings:
                # Recalculate original_scores
                scores = calculate_dot_plot(seq1, seq2, window_size, matrix_type, matrix)
                st.session_state.original_scores = scores
                st.session_state.current_settings = new_settings
            else:
                # Use existing original_scores
                scores = st.session_state.original_scores

            # Update threshold in session state
            st.session_state.threshold = threshold

            if scores is not None:
                # Apply threshold to get filtered_scores
                filtered_scores = apply_threshold(scores, threshold)
                st.session_state.filtered_scores = filtered_scores

                # Create and display heatmap
                fig = create_plot(filtered_scores, window_size, threshold, color_scale)
                st.plotly_chart(fig, use_container_width=True)

                # Add download button for the filtered scores matrix
                scores_df = pd.DataFrame(
                    filtered_scores,
                    index=[f"Seq1_{i+1}" for i in range(filtered_scores.shape[0])],
                    columns=[f"Seq2_{j+1}" for j in range(filtered_scores.shape[1])]
                )
                csv = scores_df.to_csv().encode('utf-8')
                st.download_button(
                    label="💾 Download Filtered Scores Matrix (CSV)",
                    data=csv,
                    file_name="dot_plot_filtered_scores.csv",
                    mime="text/csv",
                )

    # Display sequence information if available
    if st.session_state.sequences['seq1'] and st.session_state.sequences['seq2']:
        st.sidebar.subheader("📄 Sequence Information")
        st.sidebar.markdown(f"**Sequence 1:** {len(st.session_state.sequences['seq1'])} residues")
        st.sidebar.markdown(f"**Sequence 2:** {len(st.session_state.sequences['seq2'])} residues")

if __name__ == "__main__":
    main()
