# HCV Pipeline Deployment

This repository contains a complete deployment package for the HCV Pipeline web application.

## Contents
- Web frontend (templates and static files)
- Pipeline scripts
- Deployment script
- Requirements file

## Deployment Instructions

1. Clone this repository:
   ```bash
   git clone https://<repository-url>.git
   cd hcv_pipeline_deployment
   ```

2. Run the deployment script:
   ```bash
   ./deploy.sh
   ```

3. The script will:
   - Create a virtual environment
   - Install dependencies
   - Create a start script

4. To start the server:
   ```bash
   ./start_server.sh
   ```

5. Access the web interface at:
   ```
   http://<server-ip>:8080
   ```

## Requirements
- Python 3.6+
- Git
