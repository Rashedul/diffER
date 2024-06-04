#!/bin/bash

echo "Cleaning system to free up space..."

# Clean package cache
echo "Cleaning package cache..."
sudo apt-get clean

# Clean /tmp directory
echo "Cleaning /tmp directory..."
sudo rm -rf /tmp/*

# Clean user cache
echo "Cleaning user cache..."
rm -rf ~/.cache/*

# Remove old kernels (for Debian-based systems)
echo "Removing old kernels..."
sudo apt-get autoremove --purge

# Check disk usage
echo "Disk usage after cleanup:"
df -h

