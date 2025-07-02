document.addEventListener('DOMContentLoaded', function() {
    const pipelineForm = document.getElementById('pipeline-form');
    const runPipelineBtn = document.getElementById('run-pipeline-btn');
    const pipelineOutput = document.getElementById('pipeline-output');
    const reportsTabBtn = document.getElementById('reports-tab');
    const textReportsList = document.getElementById('text-reports-list');
    const htmlReportsList = document.getElementById('html-reports-list');
    const textViewModal = new bootstrap.Modal(document.getElementById('textViewModal'));
    const textViewContent = document.getElementById('textViewContent');
    
    // Sorting state
    let currentTextSort = 'name'; // Default to name sort
    let currentHtmlSort = 'name'; // Default to name sort
    let isTextSortAsc = true; // Default to ascending
    let isHtmlSortAsc = true; // Default to ascending

    // Initialize sort arrows
    function initSortArrows() {
        document.querySelectorAll('.sort-arrow').forEach(arrow => {
            const sortType = arrow.dataset.sort;
            const fileType = arrow.closest('.sort-btn').dataset.type;
            
            if (fileType === 'text' && sortType === currentTextSort) {
                arrow.innerHTML = isTextSortAsc ? '↑' : '↓';
            } else if (fileType === 'html' && sortType === currentHtmlSort) {
                arrow.innerHTML = isHtmlSortAsc ? '↑' : '↓';
            } else {
                arrow.innerHTML = '';
            }
        });
    }

    const textReportsPagination = document.getElementById('text-reports-pagination');
    const htmlReportsPagination = document.getElementById('html-reports-pagination');

    let allTextFiles = [];
    let allHtmlFiles = [];
    const itemsPerPage = 10; // Number of files per page

    // Function to render files for a given page
    function renderFiles(fileList, listElement, paginationElement, currentPage, fileType) {
        listElement.innerHTML = ''; // Clear current list
        paginationElement.innerHTML = ''; // Clear current pagination

        const startIndex = (currentPage - 1) * itemsPerPage;
        const endIndex = startIndex + itemsPerPage;
        const filesToRender = fileList.slice(startIndex, endIndex);

        filesToRender.forEach(file => {
            const li = document.createElement('li');
            li.className = 'list-group-item d-flex justify-content-between align-items-center';
            li.innerHTML = `
                <span>${file}</span>
                <div>
                    <button class="btn btn-sm btn-info view-file-btn" data-filename="${file}" data-filetype="${fileType}">View</button>
                    <a href="/download_file/${file}" class="btn btn-sm btn-success download-file-btn" download>Download</a>
                </div>
            `;
            listElement.appendChild(li);
        });

        // Add event listeners for view buttons
        document.querySelectorAll('.view-file-btn').forEach(button => {
            button.addEventListener('click', function() {
                const filename = this.dataset.filename;
                const filetype = this.dataset.filetype;
                if (filetype === 'text') {
                    fetch(`/view_file/${filename}`)
                        .then(response => response.text())
                        .then(content => {
                            textViewContent.textContent = content;
                            textViewModal.show();
                        })
                        .catch(error => {
                            console.error('Error viewing text file:', error);
                            alert('Error viewing file. Check console for details.');
                        });
                } else if (filetype === 'html') {
                    // Open HTML files in a new tab
                    window.open(`/view_file/${filename}`, '_blank');
                }
            });
        });

        // Render pagination controls
        const totalPages = Math.ceil(fileList.length / itemsPerPage);
        if (totalPages > 1) {
            // First button
            const firstLi = document.createElement('li');
            firstLi.className = `page-item ${currentPage === 1 ? 'disabled' : ''}`;
            firstLi.innerHTML = `<a class="page-link" href="#" data-page="1">First</a>`;
            paginationElement.appendChild(firstLi);

            // Previous button
            const prevLi = document.createElement('li');
            prevLi.className = `page-item ${currentPage === 1 ? 'disabled' : ''}`;
            prevLi.innerHTML = `<a class="page-link" href="#" data-page="${currentPage - 1}">Previous</a>`;
            paginationElement.appendChild(prevLi);

            // Page numbers
            for (let i = 1; i <= totalPages; i++) {
                const pageLi = document.createElement('li');
                pageLi.className = `page-item ${currentPage === i ? 'active' : ''}`;
                pageLi.innerHTML = `<a class="page-link" href="#" data-page="${i}">${i}</a>`;
                paginationElement.appendChild(pageLi);
            }

            // Next button
            const nextLi = document.createElement('li');
            nextLi.className = `page-item ${currentPage === totalPages ? 'disabled' : ''}`;
            nextLi.innerHTML = `<a class="page-link" href="#" data-page="${currentPage + 1}">Next</a>`;
            paginationElement.appendChild(nextLi);

            // Last button
            const lastLi = document.createElement('li');
            lastLi.className = `page-item ${currentPage === totalPages ? 'disabled' : ''}`;
            lastLi.innerHTML = `<a class="page-link" href="#" data-page="${totalPages}">Last</a>`;
            paginationElement.appendChild(lastLi);

            // Add event listeners for pagination buttons
            paginationElement.querySelectorAll('.page-link').forEach(link => {
                link.addEventListener('click', function(e) {
                    e.preventDefault();
                    const newPage = parseInt(this.dataset.page);
                    if (newPage >= 1 && newPage <= totalPages) {
                        if (listElement.id === 'text-reports-list') {
                            renderFiles(allTextFiles, textReportsList, textReportsPagination, newPage, 'text');
                        } else if (listElement.id === 'html-reports-list') {
                            renderFiles(allHtmlFiles, htmlReportsList, htmlReportsPagination, newPage, 'html');
                        }
                    }
                });
            });
        }
    }

    // Sorting functions
    function sortByName(a, b) {
        return a.localeCompare(b);
    }


    // Function to sort files based on current sort type
    function sortFiles(files, sortType, isAsc) {
        const sortedFiles = [...files];
        if (sortType === 'name') {
            sortedFiles.sort(sortByName);
        }
        return isAsc ? sortedFiles : sortedFiles.reverse();
    }

    // Add event listeners for sort buttons
    document.querySelectorAll('.sort-btn').forEach(button => {
        button.addEventListener('click', function() {
            const sortType = this.dataset.sort;
            const fileType = this.dataset.type;
            
            if (fileType === 'text') {
                if (currentTextSort === sortType) {
                    isTextSortAsc = !isTextSortAsc;
                } else {
                    currentTextSort = sortType;
                    isTextSortAsc = true;
                }
                allTextFiles = sortFiles(allTextFiles, currentTextSort, isTextSortAsc);
                renderFiles(allTextFiles, textReportsList, textReportsPagination, 1, 'text');
            } else if (fileType === 'html') {
                if (currentHtmlSort === sortType) {
                    isHtmlSortAsc = !isHtmlSortAsc;
                } else {
                    currentHtmlSort = sortType;
                    isHtmlSortAsc = true;
                }
                allHtmlFiles = sortFiles(allHtmlFiles, currentHtmlSort, isHtmlSortAsc);
                renderFiles(allHtmlFiles, htmlReportsList, htmlReportsPagination, 1, 'html');
            }
            initSortArrows();
        });
    });

    // Function to fetch and display reports
    function fetchReports() {
        fetch('/list_files')
            .then(response => response.json())
            .then(data => {
                allTextFiles = sortFiles(data.text_files, currentTextSort, isTextSortAsc);
                allHtmlFiles = sortFiles(data.html_files, currentHtmlSort, isHtmlSortAsc);

                renderFiles(allTextFiles, textReportsList, textReportsPagination, 1, 'text');
                renderFiles(allHtmlFiles, htmlReportsList, htmlReportsPagination, 1, 'html');
            })
            .catch(error => {
                console.error('Error fetching reports:', error);
                textReportsList.innerHTML = '<li class="list-group-item text-danger">Error loading reports.</li>';
                htmlReportsList.innerHTML = '<li class="list-group-item text-danger">Error loading reports.</li>';
            });
    }

    // Event listener for Run Pipeline form submission
    pipelineForm.addEventListener('submit', function(e) {
        e.preventDefault(); // Prevent default form submission

        runPipelineBtn.disabled = true; // Disable button during run
        pipelineOutput.textContent = ''; // Clear previous output

        const formData = new FormData(pipelineForm);
        const data = {};
        formData.forEach((value, key) => {
            if (key === 'keep_tmp' || key === 'keep_unmerged') {
                data[key] = formData.get(key) === 'on'; // Checkbox values are 'on' or null
            } else if (key === 'kmer_overlap_threshold' || key === 'snp_distance_threshold') {
                data[key] = parseFloat(value);
            } else {
                data[key] = value;
            }
        });

        fetch('/run_pipeline', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(data)
        })
        .then(response => {
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            const reader = response.body.getReader();
            const decoder = new TextDecoder('utf-8');
            let buffer = '';

            function readStream() {
                reader.read().then(({ done, value }) => {
                    if (done) {
                        console.log('Stream complete');
                        runPipelineBtn.disabled = false; // Re-enable button
                        return;
                    }
                    buffer += decoder.decode(value, { stream: true });
                    const lines = buffer.split('\n');
                    buffer = lines.pop(); // Keep incomplete line in buffer

                    lines.forEach(line => {
                        if (line.startsWith('data: ')) {
                            const content = line.substring(6);
                            pipelineOutput.textContent += content + '\n';
                            pipelineOutput.scrollTop = pipelineOutput.scrollHeight; // Auto-scroll
                        }
                    });
                    readStream(); // Continue reading
                }).catch(error => {
                    console.error('Stream reading error:', error);
                    pipelineOutput.textContent += `\nError: ${error.message}\n`;
                    runPipelineBtn.disabled = false; // Re-enable button
                });
            }
            readStream();
        })
        .catch(error => {
            console.error('Error starting pipeline:', error);
            pipelineOutput.textContent += `\nError: ${error.message}\n`;
            runPipelineBtn.disabled = false; // Re-enable button
        });
    });

    // Initial fetch of reports when the page loads
    fetchReports();
    initSortArrows();

    // Event listener for Reports tab activation (only re-render, don't re-fetch)
    reportsTabBtn.addEventListener('shown.bs.tab', function (event) {
        // Re-render files based on current sort state, no need to re-fetch
        renderFiles(allTextFiles, textReportsList, textReportsPagination, 1, 'text');
        renderFiles(allHtmlFiles, htmlReportsList, htmlReportsPagination, 1, 'html');
        initSortArrows(); // Ensure arrows are correct if tab was hidden and shown
    });
});
