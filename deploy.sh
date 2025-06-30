#!/bin/bash

# Deployment script for HCV Pipeline Web App

# Create deployment directory
DEPLOY_DIR="hcv_pipeline_deploy_$(date +%Y%m%d_%H%M%S)"
mkdir -p $DEPLOY_DIR

# Copy necessary files
cp -r web_frontend $DEPLOY_DIR/
cp HCV_combine_reads_V0_1.py $DEPLOY_DIR/
cp HCV_transmission_test_V0_2.py $DEPLOY_DIR/
cp HCV_Master_script_V0_1.py $DEPLOY_DIR/
cp requirements.txt $DEPLOY_DIR/

# Create virtual environment
cd $DEPLOY_DIR
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Create start script
cat > start_server.sh << 'EOL'
#!/bin/bash
source venv/bin/activate
waitress-serve --call 'web_frontend.app:create_app'
EOL

chmod +x start_server.sh

# Create README
cat > README.txt << 'EOL'
HCV Pipeline Web App Deployment

To start the server:
1. Open terminal in this directory
2. Run: ./start_server.sh

The server will run on port 8080 by default.

To access the web interface:
http://<server-ip>:8080

EOL

echo "Deployment created in directory: $DEPLOY_DIR"
echo "To deploy on another server:"
echo "1. Copy the $DEPLOY_DIR directory"
echo "2. Run ./start_server.sh"
