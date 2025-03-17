# Use an official Python image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the script into the container
COPY glycopeptide_sequence_finder_cmd.py .

# Set the default entrypoint to allow passing arguments
ENTRYPOINT ["python", "glycopeptide_sequence_finder_cmd.py"]
