#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#define MAX_BLOCKS 100000
#define MAX_SEGMENTS 10000

int active_segment_count = 0;
int next_seg_id = 0; // start at 1 to avoid confusion with default 0

// Shared range metadata for tracking views
struct shared_range_node {
    int from_seg_id;                 // ID of segment that created the view
    size_t block_start;             // start in block
    size_t block_end;               // end in block (exclusive)
    struct shared_range_node* next;
};

struct data_block {
    int16_t* samples;
    size_t length;
    struct shared_range_node* shared_views;  // linked list of shared ranges
};

struct seg_entry {
    int block_index;
    size_t block_start;
    size_t block_end;
    size_t seg_start;
    size_t seg_end;
    bool is_parent;
    size_t num_samples;
};

struct sound_seg {
    int id; // unique segment identifier
    struct seg_entry* entries;
    int num_entries;
    int capacity;
    bool has_children;
    size_t length;
};

struct data_block shared_pool[MAX_BLOCKS];
bool used_blocks[MAX_BLOCKS] = {false};

//-------------------------------------------------------------------------------------------------

int create_block(size_t initial_size) {
    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (!used_blocks[i]) {
            shared_pool[i].samples = malloc(initial_size * sizeof(int16_t));
            shared_pool[i].length = initial_size;
            shared_pool[i].shared_views = NULL;
            used_blocks[i] = true;
            return i;
        }
    }
    return -1;
}

void destroy_block_if_unused(int block_index) {
    if (!used_blocks[block_index]) return;
    if (shared_pool[block_index].shared_views != NULL) return; // still in use

    free(shared_pool[block_index].samples);
    shared_pool[block_index].samples = NULL;
    shared_pool[block_index].length = 0;
    used_blocks[block_index] = false;
}

void refresh_seg_positions(struct sound_seg* seg) {
    size_t curr = 0;
    for (int i = 0; i < seg->num_entries; i++) {
        struct seg_entry* e = &seg->entries[i];
        e->seg_start = curr;
        e->seg_end = curr + (e->block_end - e->block_start);
        e->num_samples = e->seg_end - e->seg_start;
        curr = e->seg_end;
    }
    seg->length = curr;
}

bool is_last_entry(struct sound_seg* seg, int index) {
    return index == seg->num_entries - 1;
}

void insert_entry_after(struct sound_seg* seg, int index, struct seg_entry e) {
    if (seg->num_entries >= seg->capacity) {
        seg->capacity *= 2;
        seg->entries = realloc(seg->entries, seg->capacity * sizeof(struct seg_entry));
    }
    for (int i = seg->num_entries; i > index + 1; i--) {
        seg->entries[i] = seg->entries[i - 1];
    }
    seg->entries[index + 1] = e;
    seg->num_entries++;
}

void insert_entry_at(struct sound_seg* seg, int index, struct seg_entry e) {
    if (seg->num_entries >= seg->capacity) {
        seg->capacity *= 2;
        seg->entries = realloc(seg->entries, seg->capacity * sizeof(struct seg_entry));
    }
    for (int i = seg->num_entries; i > index; i--) {
        seg->entries[i] = seg->entries[i - 1];
    }
    seg->entries[index] = e;
    seg->num_entries++;
}

void remove_entry_at(struct sound_seg* seg, int index) {
    for (int i = index; i < seg->num_entries - 1; i++) {
        seg->entries[i] = seg->entries[i + 1];
    }
    seg->num_entries--;
}

void split_entry_for_deletion(struct sound_seg* seg, int index, size_t del_start, size_t del_end) {
    struct seg_entry* orig = &seg->entries[index];

    size_t block_off1 = del_start - orig->seg_start;
    size_t block_off2 = del_end - orig->seg_start;

    struct seg_entry right = *orig;
    right.block_start = orig->block_start + block_off2;
    right.block_end = orig->block_end;
    right.seg_start = del_end;
    right.seg_end = orig->seg_end;

    orig->block_end = orig->block_start + block_off1;
    orig->seg_end = del_start;

    insert_entry_after(seg, index, right);
}

void split_entry_for_insert(struct sound_seg* seg, int index, size_t insert_pos, size_t insert_len) {
    struct seg_entry* e = &seg->entries[index];
    size_t split_offset = insert_pos - e->seg_start;
    size_t block_offset = e->block_start + split_offset;

    struct seg_entry right = *e;
    right.seg_start = insert_pos + insert_len;
    right.seg_end = e->seg_end + insert_len;
    right.block_start = block_offset;

    e->seg_end = insert_pos;
    e->block_end = block_offset;

    insert_entry_after(seg, index, right);
}

void cleanup_pool() {
    for (int i = 0; i < MAX_BLOCKS; i++) {
        if (!used_blocks[i]) continue;

        // Free remaining shared range nodes
        struct shared_range_node* node = shared_pool[i].shared_views;
        while (node) {
            struct shared_range_node* next = node->next;
            free(node);
            node = next;
        }

        shared_pool[i].shared_views = NULL;

        if (shared_pool[i].samples) {
            free(shared_pool[i].samples);
            shared_pool[i].samples = NULL;
        }

        shared_pool[i].length = 0;
        used_blocks[i] = false;
    }
}

bool is_shared_overlap(int block_index, size_t block_offset, size_t len, int from_seg_id) {
    size_t write_end = block_offset + len;
    struct shared_range_node* node = shared_pool[block_index].shared_views;

    while (node) {
        if (node->from_seg_id != from_seg_id) {
            if (!(write_end <= node->block_start || block_offset >= node->block_end)) {
                // There is an overlap with another segment's view
                return true;
            }
        }
        node = node->next;
    }

    return false; // no overlap with other segment views
}

void propagate_write_to_shared_views(int block_index, size_t write_block_start, int16_t* src, size_t src_offset, size_t len) {
    struct shared_range_node* node = shared_pool[block_index].shared_views;

    while (node) {
        size_t node_start = node->block_start;
        size_t node_end = node->block_end;

        // Overlap between [write_start, write_end) and [node_start, node_end)
        size_t overlap_start = (write_block_start > node_start) ? write_block_start : node_start;
        size_t overlap_end = ((write_block_start + len) < node_end) ? (write_block_start + len) : node_end;

        if (overlap_start < overlap_end) {
            size_t copy_len = overlap_end - overlap_start;
            size_t src_offset_local = overlap_start - write_block_start;

            memcpy(&shared_pool[block_index].samples[overlap_start],
                   &src[src_offset + src_offset_local],
                   copy_len * sizeof(int16_t));
        }

        node = node->next;
    }
}
//-------------------------------------------------------------------------------------------------
void wav_load(const char* fname, int16_t* dest){
    
    // This section declares variables for the WAV header
    FILE *file;
    char ck_ids[4]; // stores strings to check for chunk IDs (RIFF, WAVE etc.)
    int32_t filesize; 
    int32_t format_length; // chunk size 16
    int16_t format_type; // WAVE_FORMAT_PCM
    int16_t num_channels; // 1 = Mono
    int32_t sample_rate; // 8000Hz sample rate
    int32_t bytes_per_second; // sample_rate * num_channels * bits_per_sample/8
    int16_t block_align; // bits_per_sample/8 * num_channels
    int16_t bits_per_sample; // 16 bits per sample
    int32_t data_size; 


    // Declares file to read from
    file = fopen(fname, "rb");

    // Checks if file exists
    if (file == NULL) {
        printf("Failed to open file: %s", fname);
        return;
    }   


    // This section reads the WAV header from the file
    fread(ck_ids, 1, 4, file); // reads for CHUNK ID 'RIFF'
    if (ck_ids[0] != 'R' || ck_ids[1] != 'I' || ck_ids[2] != 'F' || ck_ids[3] != 'F') {
        printf("Field should read RIFF: reads %4s", ck_ids);
        fclose(file);
        return;
    }

    fread(&filesize, 4, 1, file); // reads overall file size

    fread(ck_ids, 1, 4, file); // reads for WAVE ID 'WAVE'
    if (ck_ids[0] != 'W' || ck_ids[1] != 'A' || ck_ids[2] != 'V' || ck_ids[3] != 'E') {
        printf("Field should read WAVE: reads %4s", ck_ids);
        fclose(file);
        return;
    }

    fread(ck_ids, 1, 4, file); // reads for CHUNK ID 'fmt'
    if (ck_ids[0] != 'f' || ck_ids[1] != 'm' || ck_ids[2] != 't' || ck_ids[3] != ' ') {
        printf("Field should read fmt: reads %4s", ck_ids);
        fclose(file);
        return;
    }

    fread(&format_length, 4, 1, file); // reads for format length e.g. 16 for PCM
    fread(&format_type, 2, 1, file); // reads for format type code e.g. 1
    if (format_type != 0x0001) {
        printf("Format type should be 1: reads %d", format_type);
        fclose(file);
        return;
    }

    fread(&num_channels, 2, 1, file); // reads for number of channels, 1 for Mono
    if (num_channels != 1) {
        printf("Channel should be set to 1 for Mono: reads %d", num_channels);
        fclose(file);
        return;
    }

    fread(&sample_rate, 4, 1, file); // reads for sample rate 8000Hz
    if (sample_rate != 8000) {
        printf("Sample rate should be 8000Hz: reads %d", sample_rate);
        fclose(file);
        return;
    }

    fread(&bytes_per_second, 4, 1, file);
    fread(&block_align, 2, 1, file);
    fread(&bits_per_sample, 2, 1, file);
    if (bits_per_sample != 16) {
        printf("Bits p/sample should be 16: reads %d", bits_per_sample);
        fclose(file);
        return;
    }

    fread(ck_ids, 4, 1, file);
    if (ck_ids[0] != 'd' || ck_ids[1] != 'a' || ck_ids[2] != 't' || ck_ids[3] != 'a') {
        printf("Field should read data: reads %s", ck_ids);
        fclose(file);
        return;
    }

    fread(&data_size, 4, 1, file);

    // Reads file audio data into the dest buffer
    fread(dest, 1, data_size, file);

    fclose(file);
}

// Create/write a WAV file from buffer
void wav_save(const char* fname, int16_t* src, size_t len){
    
    // This section declares variables for all the parts of the WAV header
    FILE *file;
    int32_t format_length = 16;
    int16_t format_type = 1;
    const int32_t sample_rate = 8000;
    const int16_t num_channels = 1;
    const int16_t bits_per_sample = 16;

    int32_t bytes_per_second = 16000;
    int16_t block_align = 2;
    int32_t data_size = (len*16)/8;
    int32_t filesize = 4 + (8+16) + (8+data_size);

    // Creates a file variable to write to
    file = fopen(fname, "wb");
    
    // Checks if file variable exists
    if (file == NULL) {
        printf("Failed to open file: %s\n", fname);
        return;
    }


    // This section writes the WAV header into the file
    fwrite("RIFF", 1, 4, file);
    fwrite(&filesize, 4, 1, file);
    fwrite("WAVE", 1, 4, file);
    fwrite("fmt ", 1, 4, file);

    fwrite(&format_length, 4, 1, file);
    fwrite(&format_type, 2, 1, file);
    fwrite(&num_channels, 2, 1, file);
    fwrite(&sample_rate, 4, 1, file);
    fwrite(&bytes_per_second, 4, 1, file);
    fwrite(&block_align, 2, 1, file);
    fwrite(&bits_per_sample, 2, 1, file);

    fwrite("data", 1, 4, file);
    fwrite(&data_size, 4, 1, file);


    // Writes data from src buffer into the file
    fwrite(src, 2, len, file);

    fclose(file);


    return;
}

// Initialize a new sound_seg object
struct sound_seg* tr_init() {

    struct sound_seg* seg = malloc(sizeof(struct sound_seg));
    if (!seg) return NULL;

    seg->entries = malloc(16 * sizeof(struct seg_entry)); // initial capacity
    if (!seg->entries) {
        free(seg);
        return NULL;
    }

    seg->id = next_seg_id++;  // ← Unique ID assigned here
    seg->num_entries = 0;
    seg->capacity = 16;
    seg->length = 0;
    seg->has_children = false;
    active_segment_count++;
    return seg;
}

// Destroy a sound_seg object and free all allocated memory
void tr_destroy(struct sound_seg* obj) {

    if (!obj) return;

    for (int i = 0; i < obj->num_entries; i++) {
        struct seg_entry* e = &obj->entries[i];
        if (e) {
            // Detach matching shared_range_nodes from block
            struct shared_range_node** cur = &shared_pool[e->block_index].shared_views;
            while (*cur) {
                if ((*cur)->from_seg_id == obj->id &&
                    (*cur)->block_start == e->block_start &&
                    (*cur)->block_end == e->block_end) {
                    struct shared_range_node* to_free = *cur;
                    *cur = (*cur)->next;
                    free(to_free);
                    break; // only one node per seg_entry
                } else {
                    cur = &(*cur)->next; 
                }
            }
        }
    }

    free(obj->entries);
    free(obj);

    active_segment_count--;
    if (active_segment_count == 0) {
        cleanup_pool();
    }
}

// Return the length of the segment
size_t tr_length(struct sound_seg* seg) {
    size_t total = 0;
    for (int i = 0; i < seg->num_entries; ++i) {
        total += seg->entries[i].num_samples;
    }
    return total;
}

// Read len elements from position pos into dest
void tr_read(struct sound_seg* track, int16_t* dest, size_t pos, size_t len) {
    if (!track || !dest || len == 0) return;

    size_t dest_offset = 0;
    size_t read_end = pos + len;

    for (int i = 0; i < track->num_entries && dest_offset < len; i++) {
        struct seg_entry* e = &track->entries[i];

        // Skip entries before read range
        if (e->seg_end <= pos) continue;

        // Stop if we've read past the requested range
        if (e->seg_start >= read_end) break;

        // Compute actual overlap
        size_t seg_read_start = (e->seg_start < pos) ? pos : e->seg_start;
        size_t seg_read_end   = (e->seg_end > read_end) ? read_end : e->seg_end;
        size_t to_copy = seg_read_end - seg_read_start;

        size_t block_offset = e->block_start + (seg_read_start - e->seg_start);

        memcpy(&dest[dest_offset],
               &shared_pool[e->block_index].samples[block_offset],
               to_copy * sizeof(int16_t));

        dest_offset += to_copy;
    }

}

// Write len elements from src into position pos
void tr_write(struct sound_seg* track, int16_t* src, size_t pos, size_t len) {
    if (!track || !src || len == 0) return;

    size_t write_start = pos;
    size_t write_end = pos + len;
    size_t src_offset = 0;

    // Step 1: Overwrite any overlapping entries in-place
    for (int i = 0; i < track->num_entries && src_offset < len; i++) {
        struct seg_entry* e = &track->entries[i];

        if (e->seg_end <= write_start || e->seg_start >= write_end)
            continue;

        // Compute overlap between write and entry
        size_t overlap_start = (e->seg_start > write_start) ? e->seg_start : write_start;
        size_t overlap_end   = (e->seg_end < write_end) ? e->seg_end : write_end;
        size_t overlap_len   = overlap_end - overlap_start;

        size_t block_offset  = e->block_start + (overlap_start - e->seg_start);
        size_t src_index     = overlap_start - write_start;

        memcpy(&shared_pool[e->block_index].samples[block_offset],
               &src[src_index],
               overlap_len * sizeof(int16_t));
    }

    // Step 2: Write any leftover part that didn't overlap existing entries
    // This can happen if we're extending past the end of the current segment
    size_t actual_len = tr_length(track);
    if (write_end > actual_len) {
        size_t new_len = write_end - actual_len;
        int block_index = create_block(new_len);
        if (block_index == -1) return;

        memcpy(shared_pool[block_index].samples,
               &src[len - new_len],
               new_len * sizeof(int16_t));
        shared_pool[block_index].length = new_len;

        struct seg_entry e = {
            .block_index = block_index,
            .block_start = 0,
            .block_end = new_len,
            .seg_start = 0, // filled by refresh
            .seg_end = 0,
            .is_parent = false,
            .num_samples = new_len
        };

        insert_entry_at(track, track->num_entries, e);
        refresh_seg_positions(track);
    }

}

// Delete a range of elements from the track
bool tr_delete_range(struct sound_seg* track, size_t pos, size_t len) {
    if (!track || len == 0 || track->num_entries == 0) return true;

    size_t del_start = pos;
    size_t del_end = pos + len;

    // REQ 3.1: Deletion must respect shared parent content
    for (int i = 0; i < track->num_entries; ++i) {
        struct seg_entry* e = &track->entries[i];
        if (e->is_parent) {
            if (del_start < e->seg_end && del_end > e->seg_start) {
                return false;
            }
        }
    }

    // Traverse and modify entries
    int i = 0;
    while (i < track->num_entries) {
        struct seg_entry* e = &track->entries[i];

        if (e->seg_end <= del_start) {
            i++; continue;
        }

        if (e->seg_start >= del_end) {
            break;
        }

        size_t overlap_start = (del_start > e->seg_start) ? del_start : e->seg_start;
        size_t overlap_end   = (del_end < e->seg_end) ? del_end : e->seg_end;

        if (overlap_start > e->seg_start && overlap_end < e->seg_end) {
            // Split case — convert to left and right
            struct seg_entry right = *e;

            right.seg_start = overlap_end;
            right.block_start += (overlap_end - e->seg_start);
            right.seg_end = e->seg_end;
            right.block_end = e->block_end;
            right.num_samples = right.seg_end - right.seg_start;

            e->seg_end = overlap_start;
            e->block_end = e->block_start + (overlap_start - e->seg_start);
            e->num_samples = e->seg_end - e->seg_start;

            insert_entry_after(track, i, right);
            i += 2;
            continue;
        }

        if (overlap_start == e->seg_start && overlap_end < e->seg_end) {
            // Trim left
            size_t delta = overlap_end - e->seg_start;
            e->block_start += delta;
            e->seg_start = overlap_end;
            e->num_samples = e->seg_end - e->seg_start;
            i++;
            continue;
        }

        if (overlap_start > e->seg_start && overlap_end == e->seg_end) {
            // Trim right
            e->seg_end = overlap_start;
            e->block_end = e->block_start + (overlap_start - e->seg_start);
            e->num_samples = e->seg_end - e->seg_start;
            i++;
            continue;
        }

        if (overlap_start <= e->seg_start && overlap_end >= e->seg_end) {
            // Full delete
            remove_entry_at(track, i);
            continue;
        }

        i++;
    }

    refresh_seg_positions(track);
    return true;
}

// Returns a string containing <start>,<end> ad pairs in target
char* tr_identify(struct sound_seg* target, struct sound_seg* ad) {
    if (!target || !ad) return strdup("");

    size_t target_len = tr_length((struct sound_seg*)target);
    size_t ad_len = tr_length((struct sound_seg*)ad);

    if (ad_len == 0 || ad_len > target_len) return strdup("");

    // Allocate buffers dynamically
    int16_t* ad_buf = malloc(ad_len * sizeof(int16_t));
    int16_t* target_buf = malloc(ad_len * sizeof(int16_t));
    if (!ad_buf || !target_buf) {
        free(ad_buf);
        free(target_buf);
        return strdup("");
    }

    // Read ad segment into buffer
    tr_read((struct sound_seg*)ad, ad_buf, 0, ad_len);

    // Autocorrelation reference value
    long long ad_auto_corr = 0;
    for (size_t i = 0; i < ad_len; i++) {
        ad_auto_corr += (long long)ad_buf[i] * ad_buf[i];
    }

    double threshold = 0.95 * ad_auto_corr;

    // Result string
    char* result = malloc(1);
    result[0] = '\0';
    size_t result_len = 0;

    for (size_t start = 0; start <= target_len - ad_len; start++) {
        tr_read((struct sound_seg*)target, target_buf, start, ad_len);

        long long cross_corr = 0;
        for (size_t i = 0; i < ad_len; i++) {
            cross_corr += (long long)target_buf[i] * ad_buf[i];
        }

        if (cross_corr >= threshold) {
            size_t end = start + ad_len - 1;

            char buffer[50];
            snprintf(buffer, sizeof(buffer), "%zu,%zu\n", start, end);

            size_t buffer_len = strlen(buffer);
            char* new_result = realloc(result, result_len + buffer_len + 1);
            if (!new_result) break;

            result = new_result;
            strcpy(result + result_len, buffer);
            result_len += buffer_len;

            // skip to avoid overlap
            start = end;
        }
    }

    if (result_len > 0 && result[result_len - 1] == '\n') {
        result[result_len - 1] = '\0';
    }

    free(ad_buf);
    free(target_buf);
    return result;
}

// Insert a portion of src_track into dest_track at position destpos
void tr_insert(struct sound_seg* src_track,
    struct sound_seg* dest_track,
    size_t destpos, size_t srcpos, size_t len) {

    if (!src_track || !dest_track || len == 0) return;

    size_t src_len = tr_length(src_track);
    size_t dest_len = tr_length(dest_track);
    if (srcpos + len > src_len || destpos > dest_len) return;

    struct seg_entry* snapshot = malloc(sizeof(struct seg_entry) * src_track->num_entries);
    int snap_count = 0;
    size_t inserted = 0;

    // === STEP 1: Collect actual block slices from source ===
    for (int i = 0; i < src_track->num_entries && inserted < len; i++) {
        struct seg_entry* e = &src_track->entries[i];

        // Skip entries that don't intersect with the desired source range
        if (e->seg_end <= srcpos || e->seg_start >= srcpos + len)
            continue;

        // Compute overlap
        size_t view_start = (srcpos > e->seg_start) ? srcpos : e->seg_start;
        size_t view_end   = ((srcpos + len) < e->seg_end) ? (srcpos + len) : e->seg_end;
        size_t view_len   = view_end - view_start;

        // Compute real block offset inside the block (even if it's from a view)
        size_t block_offset = e->block_start + (view_start - e->seg_start);
        size_t block_end    = block_offset + view_len;

        struct seg_entry view = {
            .block_index = e->block_index,
            .block_start = block_offset,
            .block_end = block_end,
            .seg_start = 0,
            .seg_end = 0,
            .is_parent = false,
            .num_samples = view_len
        };

        // Mark original entry as parent
        e->is_parent = true;

        // Track view in shared range
        struct shared_range_node* node = malloc(sizeof(struct shared_range_node));
        node->from_seg_id = dest_track->id;
        node->block_start = block_offset;
        node->block_end = block_end;
        node->next = shared_pool[e->block_index].shared_views;
        shared_pool[e->block_index].shared_views = node;

        snapshot[snap_count++] = view;
        inserted += view_len;
    }

    // === STEP 2: Split target entry at insert point, if needed ===
    int insert_at = 0;
    while (insert_at < dest_track->num_entries &&
           dest_track->entries[insert_at].seg_end <= destpos) {
        insert_at++;
    }

    if (insert_at < dest_track->num_entries) {
        struct seg_entry* orig = &dest_track->entries[insert_at];
        if (destpos > orig->seg_start && destpos < orig->seg_end) {
            size_t offset = destpos - orig->seg_start;
            size_t block_split = orig->block_start + offset;

            struct seg_entry right = *orig;
            right.block_start = block_split;
            right.seg_start = 0;
            right.seg_end = 0;
            right.num_samples = orig->block_end - block_split;

            orig->block_end = block_split;
            orig->num_samples = offset;

            insert_entry_after(dest_track, insert_at, right);
            insert_at++;
        }
    }

    // === STEP 3: Shift right and insert new views ===
    for (int i = dest_track->num_entries - 1; i >= insert_at; i--) {
        dest_track->entries[i + snap_count] = dest_track->entries[i];
    }
    for (int i = 0; i < snap_count; i++) {
        dest_track->entries[insert_at + i] = snapshot[i];
    }
    dest_track->num_entries += snap_count;

    // === STEP 4: Refresh timeline ===
    refresh_seg_positions(dest_track);
    free(snapshot);
}



int main() {
    
    
}