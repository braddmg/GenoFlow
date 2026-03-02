# Convert Docker image to Singularity container

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_header() {
    echo -e "${BLUE}$1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}→ $1${NC}"
}

# Check prerequisites
check_singularity() {
    if ! command -v singularity &> /dev/null; then
        print_error "Singularity not found. Install Singularity and try again."
        echo ""
        echo "Installation instructions:"
        echo "  Ubuntu/Debian: sudo apt-get install -y singularity-container"
        echo "  CentOS/RHEL:   sudo yum install -y singularity"
        echo "  macOS (via Vagrant): https://sylabs.io/guides/latest/user-guide/"
        exit 1
    fi
    print_success "Singularity found: $(singularity --version)"
}

check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker not found. Install Docker first."
        exit 1
    fi
    print_success "Docker found: $(docker --version)"
}

check_docker_daemon() {
    if ! docker ps &> /dev/null; then
        print_error "Docker daemon not running."
        exit 1
    fi
    print_success "Docker daemon running"
}

# Build from definition file
build_from_def() {
    print_header "Building Singularity from definition file"
    
    if [ ! -f Singularity.def ]; then
        print_error "Singularity.def not found in current directory"
        exit 1
    fi
    
    output_file="${1:-genflow.sif}"
    
    print_info "Building: $output_file"
    print_info "This may take 10-30 minutes..."
    echo ""
    
    if sudo singularity build "$output_file" Singularity.def; then
        print_success "Singularity image built: $output_file"
        ls -lh "$output_file"
        return 0
    else
        print_error "Build failed. Check error messages above."
        return 1
    fi
}

# Convert Docker to Singularity
convert_docker_to_singularity() {
    print_header "Converting Docker image to Singularity"
    
    docker_image="${1:-genflow:latest}"
    output_file="${2:-genflow.sif}"
    
    print_info "Docker image: $docker_image"
    print_info "Output file: $output_file"
    echo ""
    
    # Check if Docker image exists
    if ! docker images | grep -q "$docker_image"; then
        print_error "Docker image not found: $docker_image"
        echo "Available Docker images:"
        docker images | grep genflow || echo "  (none)"
        return 1
    fi
    
    print_info "This may take several minutes..."
    echo ""
    
    # Create temporary directory
    temp_dir=$(mktemp -d)
    print_info "Using temporary directory: $temp_dir"
    
    # Export Docker image to tar
    print_info "Exporting Docker image..."
    docker save "$docker_image" -o "$temp_dir/docker.tar"
    
    # Convert Docker tar to Singularity image
    print_info "Converting to Singularity format..."
    if singularity build --force "$output_file" "docker-archive://$temp_dir/docker.tar"; then
        print_success "Singularity image created: $output_file"
        ls -lh "$output_file"
        
        # Cleanup
        rm -rf "$temp_dir"
        print_success "Temporary files cleaned up"
        return 0
    else
        print_error "Conversion failed"
        rm -rf "$temp_dir"
        return 1
    fi
}

# Test Singularity image
test_singularity() {
    sif_file="${1:-genflow.sif}"
    
    print_header "Testing Singularity image"
    
    if [ ! -f "$sif_file" ]; then
        print_error "Singularity image not found: $sif_file"
        return 1
    fi
    
    print_info "File: $sif_file ($(ls -lh "$sif_file" | awk '{print $5}'))"
    echo ""
    
    print_info "Test 1: Container shell"
    if singularity exec "$sif_file" echo "✓ Shell works"; then
        print_success "Test 1 passed"
    else
        print_error "Test 1 failed"
        return 1
    fi
    
    print_info "Test 2: Conda environments"
    if singularity exec "$sif_file" bash -c "conda env list | grep genflow" &>/dev/null; then
        print_success "Test 2 passed"
    else
        print_error "Test 2 failed"
        return 1
    fi
    
    print_info "Test 3: Snakemake"
    if singularity exec "$sif_file" snakemake --version &>/dev/null; then
        version=$(singularity exec "$sif_file" snakemake --version)
        print_success "Test 3 passed: Snakemake $version"
    else
        print_error "Test 3 failed"
        return 1
    fi
    
    print_info "Test 4: Prokka"
    if singularity exec "$sif_file" bash -c "conda run -n genflow-prokka prokka --version" &>/dev/null; then
        print_success "Test 4 passed"
    else
        print_error "Test 4 failed"
        return 1
    fi
    
    echo ""
    print_success "All tests passed!"
}

# Show usage
show_usage() {
    cat << 'EOF'
GenFlow Docker → Singularity Converter

Usage:
    ./docker-to-singularity.sh [COMMAND] [OPTIONS]

Commands:
    build-def [output.sif]      Build from Singularity.def file
    convert [image] [output]    Convert Docker image to Singularity
    test [image.sif]            Test Singularity image
    help                        Show this help message

Examples:
    # Build from definition file
    ./docker-to-singularity.sh build-def genflow.sif
    
    # Convert existing Docker image
    ./docker-to-singularity.sh convert genflow:latest genflow.sif
    
    # Test the resulting image
    ./docker-to-singularity.sh test genflow.sif

Workflow:
    1. Build Docker image:
       docker-compose build
    
    2. Convert to Singularity:
       ./docker-to-singularity.sh convert genflow:latest genflow.sif
    
    3. Test Singularity image:
       ./docker-to-singularity.sh test genflow.sif
    
    4. Use the container:
       singularity run genflow.sif snakemake --cores 8 --use-conda \
           --config fasta='Data/genome.fasta'

Notes:
    - Requires sudo for Singularity build
    - Conversion may take 10-30 minutes
    - Output image size is typically 8-15 GB
    - Best performed with 50+ GB free disk space

EOF
}

# Main script
main() {
    print_header "GenFlow Docker → Singularity Converter"
    echo ""
    
    # Parse arguments
    command="${1:-build-def}"
    
    case "$command" in
        build-def)
            check_singularity
            build_from_def "$2"
            if [ $? -eq 0 ] && [ -n "$2" ]; then
                echo ""
                read -p "Test the image now? (y/n): " test_now
                if [ "$test_now" = "y" ]; then
                    test_singularity "$2"
                fi
            fi
            ;;
        convert)
            check_docker
            check_docker_daemon
            check_singularity
            
            docker_img="${2:-genflow:latest}"
            output="${3:-genflow.sif}"
            
            convert_docker_to_singularity "$docker_img" "$output"
            if [ $? -eq 0 ]; then
                echo ""
                read -p "Test the image now? (y/n): " test_now
                if [ "$test_now" = "y" ]; then
                    test_singularity "$output"
                fi
            fi
            ;;
        test)
            sif="${2:-genflow.sif}"
            test_singularity "$sif"
            ;;
        help|--help|-h)
            show_usage
            ;;
        *)
            print_error "Unknown command: $command"
            echo ""
            show_usage
            exit 1
            ;;
    esac
}

main "$@"
