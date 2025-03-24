# Use an official Python image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Create the necessary directory structure
RUN mkdir -p /app/src /app/data/test_proteomes /app/data/digested_peptide_library /app/data/digested_glycopeptide_library /app/data/logs

# Copy the source code
COPY src/glycopeptide_sequence_finder_cmd.py /app/src/

# Set the default entrypoint to allow passing arguments
ENTRYPOINT ["python", "/app/src/glycopeptide_sequence_finder_cmd.py"]

# Add a default CMD to show usage instructions
CMD ["--help"]
